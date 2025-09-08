/*
 * TBamFile.cpp
 *
 *  Created on: May 18, 2020
 *      Author: wegmannd
 */

#include "TBamFile.h"
#include "TAlignment.h"
#include "TOutputBamFile.h"
#include "api/BamWriter.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Strings/stringConversions.h"
#include "coretools/algorithms.h"
#include "genometools/GenomePositions/TGenomePosition.h"

namespace BAM{
using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::user_assert;

namespace impl {
std::string millionReadsRead(size_t N) { return coretools::str::toStringWithPrecision((double) N / 1000000, 1); }
}

void TBamFile::curFilterOut(){
	_filters.filterOut(FilterType::External, _read.bamAlignment.Name, _read.bamAlignment.IsReverseStrand(), _read.readGroupID, refID());
}

void TBamFile::filterOut(const TAlignment & Alignment){
	_filters.filterOut(FilterType::External, Alignment.name(), Alignment.isReverseStrand(), Alignment.readGroupId(), Alignment.refID());

}

void TBamFile::setExternalFilterReason(std::string_view reason){
	_filters.enable(FilterType::External, reason);
}

void TBamFile::_fillSamHeader(){
	//Note: chromosomes and read groups are in separate objects
	_samHeader.set(_bamHeader.Version, _bamHeader.SortOrder, _bamHeader.GroupOrder, "none");

	//add programs
	for(auto it = _bamHeader.Programs.Begin(); it != _bamHeader.Programs.End(); ++it){
		_samHeader.addProgram(it->ID, it->Name, it->CommandLine, "", it->Version);
	}

	//add links among programs
	for(auto it = _bamHeader.Programs.Begin(); it != _bamHeader.Programs.End(); ++it){
		if(it->HasPreviousProgramID()){
			_samHeader.addPreviousProgramInChain(it->ID, it->PreviousProgramID);
		}
	}

	//add comments
	for(auto& c : _bamHeader.Comments){
		_samHeader.addComment(c);
	}
}

void TBamFile::_fillChromosomes(){
	user_assert(_bamHeader.Sequences.Size() > 0, "No chromosomes present in BAM header!");

	// make sure object is empty
	_chromosomes.clear();

	//copy from BamHeader
	for(BamTools::SamSequenceIterator chrIt=_bamHeader.Sequences.Begin(); chrIt!=_bamHeader.Sequences.End(); ++chrIt){
		_chromosomes.appendChromosome(chrIt->Name, coretools::str::fromString<uint64_t>(chrIt->Length));
	}
}

void TBamFile::_fillReadGroups(){
	//make sure they are empty
	_readGroups.clear();

	//now add one by one
	//TODO : not nice how it works, but implemented this way to ensure TReadGroups does not depend on bamtools
	for(auto it = _bamHeader.ReadGroups.Begin(); it != _bamHeader.ReadGroups.End(); ++it){
		//add read group
		TReadGroup& rg = _readGroups.add(it->ID);

		//now copy rest
		rg.description_DS = it->Description;
		rg.flowOrder_FO = it->FlowOrder;
		rg.keySequence_KS = it->KeySequence;
		rg.library_LB = it->Library;
		rg.platformUnit_PU = it->PlatformUnit;
		rg.predictedInsertSize_PI = it->PredictedInsertSize;
		rg.productionDate_DT = it->ProductionDate;
		rg.program_PG = it->Program;
		rg.sample_SM = it->Sample;
		rg.sequencingCenter_CN = it->SequencingCenter;
		rg.sequencingTechnology_PL = it->SequencingTechnology;
	}
}

//--------------------------------------------------------
// Functions for reading
//--------------------------------------------------------
TBamFile::TBamFile(std::string_view Filename, size_t ID, bool EnableFilters) : _filename(Filename), _ID(ID), _filters(EnableFilters){

	//open BAM file
	logfile().list("Opening BAM file '", _filename, "'.");
	if (!_bamReader.Open(_filename))
		throw coretools::TUserError("Failed to open BAM file '", Filename, "'!");

	//load or create index file
	const std::string fnIndex1 = std::string(Filename).append(".bai");
	Filename.remove_suffix(4);
	const std::string fnIndex2 = std::string(Filename).append(".bai");
	if (std::filesystem::exists(fnIndex1)) {
		logfile().list("Opening BAM index file '", fnIndex1, "'.");
		if(!_bamReader.OpenIndex(fnIndex1))
			throw coretools::TUserError("Failed to open BAM index file '", fnIndex1, "'!");
	}
	else if (std::filesystem::exists(fnIndex2)) {
		logfile().list("Opening BAM index file '", fnIndex2, "'.");
		if (!_bamReader.OpenIndex(fnIndex2))
			throw coretools::TUserError("Failed to open BAM index file '", fnIndex2, "'!");
	} else {
		logfile().list("Creating BAM index file '", fnIndex1, "'.");
		if (!_bamReader.CreateIndex())
			throw coretools::TUserError("Failed to create BAM index file '", fnIndex1, "'!");
	}

	//initialize bam stuff
	_bamHeader = _bamReader.GetHeader();

	_fillSamHeader();

	//initialize read groups
	_fillReadGroups();
	_numNotAligned.resize(_readGroups.size());

	//initialize chromosomes
	_fillChromosomes();

	//resize alignmentCounter
	_numAlignmentReadPerReadGroupPerChromosome.resize(_readGroups.size());
	_numAlignmentReadPerReadGroupPerChromosome.resizeDistributions(_chromosomes.size());

	//get file size
	const size_t lastChromRefID = _chromosomes.size() - 1;
	int pos               = _chromosomes[lastChromRefID].length() - 1;
	BamTools::BamAlignment bamAlignment;
	do {
	_bamReader.Jump(lastChromRefID, pos);
	pos -= _step;
	} while (!_bamReader.GetNextAlignmentCore(bamAlignment) && pos > 0);
	_fileSize = _bamReader.Tell();
	_bamReader.Rewind();


	// check if we add a RG for reads without one
	if(parameters().exists("keepReadsWithoutRG")){
		std::string name = parameters().get("keepReadsWithoutRG");
		if(name.empty()){
			name = "RGForReadsWithoutReadGroup";
		}		
		logfile().list("Will put reads without read group into read group '", name, "'. (parameter 'keepReadsWithoutRG')");
		_readGroups.addReadGroupForReadsWithoutReadGroup(name);
	} else {
		logfile().list("Will filter out all reads without read group. (keep with 'keepReadsWithoutRG')");
	}	
	
	_filters.resize(_readGroups.size(), _chromosomes.size(), _filename);

	// Set Limits:
	//number of reads
	if(parameters().exists("limitReads")){
		_maxNumReadsToRead = parameters().get<uint64_t>("limitReads");
		logfile().list("Will limit the analysis to the first ", _maxNumReadsToRead, " reads in the BAM file.");
	}

	//limit chromosomes?
	_chromosomes.limitAndSetPloidy();

	//limit read groups
	if(parameters().exists("readGroup")){
		_readGroups.filterReadGroups(parameters().get<std::string>("readGroup"));
		logfile().startIndent("Will limit analysis to the following read groups:");
		_readGroups.printReadgroupsInUse();
		logfile().endIndent();
		_filters.enable(FilterType::ReadGroup, "Read group not in use");
	} else {
		_filters.disable(FilterType::ReadGroup);
		logfile().list("Will not filter Readgroups. (use 'readGroup')");
	}

	constexpr std::string_view downsample = "downsampleReads";
	_downProb = parameters().get(downsample, coretools::P(0.));
	if (_downProb > 0.) {
		logfile().list("Will downsample reads with probability ", _downProb, ".(parameter '", downsample, "')");
	} else {
		logfile().list("Will not downsample reads.(use '", downsample, "')");
	}


	constexpr std::string_view sDupReset = "resetDuplicates";
	_resetDuplicates = parameters().exists(sDupReset);
	if (_resetDuplicates) {
		logfile().list("Will reset existing duplicates tags. (parameter '", sDupReset, "')");
	} else {
		logfile().list("Will keep existing duplicates tags. (use '", sDupReset, "')");
	}

	constexpr std::string_view sDup = "markDuplicates";
	if (parameters().exists(sDup)) {
		if (parameters().get(sDup).empty()) {
			_maxDupOverlap = 0;
		} else {
			_maxDupOverlap = parameters().get<size_t>(sDup);
		}
		std::string agn = "";
		if (parameters().exists("RGagnostic")) {
			_old.resize(1);
			agn = "all readgroups";
		} else {
			_old.resize(_readGroups.size());
			agn = "same readgroup";
		}
		logfile().list("Will identify and mark duplicates on ", agn, " where start and end do not differ more than ",
					   _maxDupOverlap, ". (parameter '", sDup, "')");
	} else {
		logfile().list("Will not identify and mark duplicates. (use '", sDup, "')");
		_maxDupOverlap = _nope;
	}
}

void TBamFile::setFilters(const TBamFilters& Filters) {
	_filters.clone(Filters);
}

bool TBamFile::_readNextAlignmentFromFile(){
	using coretools::instances::randomGenerator;
	for (;;) {
		if (!_bamReader.GetNextAlignment(_read.bamAlignment)) { return false; }
		if (_downProb == 0. || randomGenerator().getRand() < _downProb) { break; }
		// Else downsaple alignment
		++_numDownsampled;
	}

	++_numAlignmentRead;

	//store current read group ID
	std::string readGroup;
	_read.bamAlignment.GetTag("RG", readGroup);
	_read.readGroupID = _readGroups.getId(readGroup);

	return true;
}

bool TBamFile::_identifyDuplicate() {
	const auto ID    = _old.size() > 1 ? _read.readGroupID : 0;
	const auto &pNow = curPosition();
	const auto &pOld = _old[ID].position;
	if (pNow.refID() != pOld.refID()) return false;

	const auto dStart = pNow.position() - pOld.position(); // sorted positions, always positive
	if (dStart > _maxDupOverlap) return false;

	const size_t dEnd =
		std::abs(int64_t(pNow.position() + _read.fragmentLength()) - int64_t(pOld.position() + _old[ID].length));
	if (dEnd > _maxDupOverlap) return false;

	_read.bamAlignment.SetIsDuplicate(true);
	return true;
}

bool TBamFile::_applyFilters() {
	// MappedLength filter is always set
	if (!_filters.pass(FilterType::MappedLength, _read.cigar.lengthMapped(), _read.bamAlignment.Name,
					   _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID())) {
		return false;
	}
	if (!_filters.enabled()) {
		return  true;
	} 


	float PMDS = -1e20;

	// apply regular filters
	return _filters.pass(FilterType::Duplicate, !_read.bamAlignment.IsDuplicate(), _read.bamAlignment.Name,
	                     _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
		   _filters.pass(FilterType::SoftClippedRation,
						 double(_read.cigar.lengthSoftClipped()) / _read.cigar.lengthRead() <= _filters.softClipRation(),
	                     _read.bamAlignment.Name, _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::PMDS, !_read.bamAlignment.GetTag("DS", PMDS) || PMDS < _filters.PMDSmax(),
	                     _read.bamAlignment.Name, _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::ImproperPairs, !_read.bamAlignment.IsPaired() || _read.bamAlignment.IsProperPair(),
	                     _read.bamAlignment.Name, _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::Unmapped, _read.bamAlignment.IsMapped(), _read.bamAlignment.Name,
	                     _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::FailedQC, !_read.bamAlignment.IsFailedQC(), _read.bamAlignment.Name,
	                     _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::Secondary, _read.bamAlignment.IsPrimaryAlignment(), _read.bamAlignment.Name,
	                     _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::Supplementary, !_read.bamAlignment.IsSupplementary(), _read.bamAlignment.Name,
	                     _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::ReadGroup, _readGroups.readGroupInUse(_read.readGroupID), _read.bamAlignment.Name,
	                     _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::FwdStrand, _read.bamAlignment.IsReverseStrand(), _read.bamAlignment.Name,
	                     _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::RevStrand, !_read.bamAlignment.IsReverseStrand(), _read.bamAlignment.Name,
	                     _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::FirstMate, _read.bamAlignment.IsFirstMate(), _read.bamAlignment.Name,
	                     _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::SecondMate, _read.bamAlignment.IsSecondMate(), _read.bamAlignment.Name,
	                     _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::MappingQuality, (size_t)_read.bamAlignment.MapQuality, _read.bamAlignment.Name,
	                     _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::Blacklist, !_filters.blacklist().isInBlacklist(_read.bamAlignment.Name),
	                     _read.bamAlignment.Name, _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::ReadLength, _read.cigar.lengthRead(), _read.bamAlignment.Name,
	                     _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::FragmentLength, _read.fragmentLength(), _read.bamAlignment.Name,
	                     _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID()) &&
	       _filters.pass(FilterType::LongerThanFragment,
	                     !_read.bamAlignment.IsProperPair() ||
	                         abs(_read.bamAlignment.InsertSize) >= static_cast<int32_t>(_read.cigar.lengthAligned()),
	                     _read.bamAlignment.Name, _read.bamAlignment.IsSecondMate(), _read.readGroupID, refID());
}

bool TBamFile::readNextAlignment(){
	if (_numAlignmentRead >= _maxNumReadsToRead) { return false; }

	//store previous position
	_previousAlignmentPosition = _read.position;


	//check if it has no read group
	bool pass = false;
	const bool justStarted = atStart();
	while (!pass) {
		pass = true;
		// get next alignment
		if (!_readNextAlignmentFromFile()) { return false; }

		if (_read.readGroupID == TReadGroups::noReadGroupId) {
			++_numNoReadGroup;
			pass = false;
		}

		// check if it is unaligned (refID < 0), in which case we read until the first aligned read
		if (_read.bamAlignment.RefID < 0) {
			++_numNotAligned[_read.readGroupID];
			pass = false;
		}
	}
	_read.position.move(_read.bamAlignment.RefID, _read.bamAlignment.Position);

	//check if BAM file is sorted
	user_assert(justStarted || _read.position >= _previousAlignmentPosition, "BAM file must be sorted by position! Alignment '",
				_read.bamAlignment.Name, "' is at position ", _read.bamAlignment.Position,
				", which is before the position of the previous alignment (", _previousAlignmentPosition.position(),
				")");

	//check if chromosome changed
	if (justStarted || chrChanged()) {
		if (refID() >= _chromosomes.size()) {
			throw coretools::TUserError("Chromosome with refID ", _read.bamAlignment.RefID,
										" is missing from BAM header!");
		}

		if(!curChromosome().inUse()) {
			do {
				_read.position.move(refID() + 1, 0);
				if(refID() >= _chromosomes.size()){
					return false;
				}
			} while (!curChromosome().inUse());
			//jump reader and read first alignment
			if (!jump(curChromosome().from())) return false;
			return readNextAlignment();
		}
	}

	//get current position, clear CIGAR and update counter
	_read.cigar.clear(); //needs to be cleared here to be empty in case of alignments that are unaligned


	// update per read group counter
	if(_read.readGroupID != TReadGroups::noReadGroupId){
		_numAlignmentReadPerReadGroupPerChromosome.add(_read.readGroupID, refID());
	}

	//parse CIGAR
	for(auto& it : _read.bamAlignment.CigarData){
		_read.cigar.add(it.Type, it.Length);
	}

	// duplicates
	if (_resetDuplicates) _read.bamAlignment.SetIsDuplicate(false);
	if (_maxDupOverlap != _nope && !_read.isDuplicate()) {
		_numIdentifiedDuplicates       += _identifyDuplicate();
		if (_old.size() > 1) {
			_old[_read.readGroupID].position = curPosition();
			_old[_read.readGroupID].length   = _read.fragmentLength();
		} else {
			_old.front().position = curPosition();
			_old.front().length   = _read.fragmentLength();
		}
	}

	//apply filters
	_read.QCFiltersPassed         = _applyFilters();
	_numAlignmentsPassedQC += _read.QCFiltersPassed;

	return true;
}

bool TBamFile::readNextAlignmentThatPassesFilters(){
	_read.QCFiltersPassed = false;
	while(!_read.QCFiltersPassed){
		if(!readNextAlignment()){
			return false;
		}
	}
	return true;
}

void TBamFile::fill(TAlignment & alignment) const{
	alignment.fill(_read.bamAlignment.Name,
				   TSamFlags(_read.bamAlignment.AlignmentFlag),
				   _read.bamAlignment.RefID,
				   _read.bamAlignment.Position,
				   _read.bamAlignment.MapQuality,
				   _read.cigar,
				   _read.bamAlignment.MateRefID,
				   _read.bamAlignment.MatePosition,
				   _read.bamAlignment.InsertSize,
				   _read.bamAlignment.QueryBases,
				   _read.bamAlignment.Qualities,
				   _ID,
				   _read.readGroupID);
}

bool TBamFile::jump(const genometools::TGenomePosition Position){
	_previousAlignmentPosition.clear();
	_read.position.clear();
	return _bamReader.Jump(Position.refID(), Position.position());
}

bool TBamFile::jump(size_t RefID){
	return jump(genometools::TGenomePosition(RefID, 0));
}

double TBamFile::averageDepth() {
	if (!atStart()) {
		logfile().warning("Calculating average depth resets bam file!");
	}
	_bamReader.Rewind();

	size_t chrLenght      = 0;
	for (const auto& chr: _chromosomes) {
		if (chr.ploidy() == 1) {
			chrLenght += chr.length()/2;
		} else {
			chrLenght += chr.length();
		}
	}

	size_t alnLength = 0;
	while (readNextAlignmentThatPassesFilters()) {
		alnLength += _read.cigar.lengthAligned();
	}

	_bamReader.Rewind();
	_read.position.clear();
	_numAlignmentRead = 0;

	return double(alnLength)/chrLenght;
}

//--------------------------------------------------------
// Functions for writing
//--------------------------------------------------------

void TBamFile::writeCurAlignment(TOutputBamFile & out){
	out.writeAlignment(_read.bamAlignment);
}

//-----------------------------------------------------
// Reporting
//-----------------------------------------------------
void TBamFile::_writeFilteringStats(std::string_view outputName) const {
	std::string filename = std::string(outputName).append("_filterSummary.txt");
	coretools::instances::logfile().listFlush("Writing general filter counts to '" + filename + "' ...");

	//creating header
	std::vector<std::string> header;
	header.push_back("ReadGroup");
	header.push_back("No_read_group");
	header.push_back("Not_aligned");
	_filters.fillHeader(header);
	coretools::TOutputFile out(filename, header, "\t");

	out.write("allReadGroups", _numNoReadGroup, coretools::containerSum(_numNotAligned));
	_filters.printCombinedCounts(out);
	out.endln();

	//writes numbers of removed reads for each applied filter per read group, also lists filters if no reads were removed
	for (size_t rg = 0; rg < _readGroups.size(); rg++){
		out.write(_readGroups.getName(rg), 0, rg < _numNotAligned.size() ? _numNotAligned[rg]: 0);
		_filters.printCounts(out, rg);
		out.endln();
	}
	out.close();
	coretools::instances::logfile().done();
}

void TBamFile::printSummary(std::string_view outputName) const {
	logfile().startIndent("Summary of parsed reads from BAM file '" + _filename + "':");
	logfile().list("Total number of reads read: ", _numAlignmentRead);
	logfile().list("Reads without read group: ", _numNoReadGroup, " (", coretools::str::toPercentString(_numNoReadGroup, _numAlignmentRead, 3), "%)");
	logfile().list("Reads identified and marked as  duplicates: ", _numIdentifiedDuplicates, " (", coretools::str::toPercentString(_numIdentifiedDuplicates, _numAlignmentRead, 3), "%)");
	logfile().list("Reads that passed filters: ", _numAlignmentsPassedQC, " (", coretools::str::toPercentString(_numAlignmentsPassedQC, _numAlignmentRead, 3), "%)");
	const auto numFiltered = _numAlignmentRead - _numAlignmentsPassedQC;
	logfile().list("Reads that were filtered out: ", numFiltered, " (" + coretools::str::toPercentString(numFiltered, _numAlignmentRead, 3), "%)");

	//write counts of filtered reads for each read group to _filterSummary.txt file
	_writeFilteringStats(outputName);

	//print counts of filtered reads for each read group to terminal, doesn't list filters if no reads were removed
	for (size_t rg = 0; rg < _readGroups.size(); rg++){
		//logfile().newLine();
		logfile().list("Number of reads filtered from read group '" + coretools::str::toString(_readGroups.getName(rg)) + "':");
		logfile().addIndent();
		bool summaryWritten = false;
		if (rg < _numNotAligned.size() && _numNotAligned[rg] > 0){
			logfile().list("Not aligned: ", _numNotAligned[rg]);
			summaryWritten = true;
		} 
		summaryWritten = summaryWritten || _filters.summary(numFiltered, rg);
		if(!summaryWritten){
			logfile().list("All reads passed filters");
		}
		logfile().endIndent();
	}

	logfile().endIndent();
	logfile().endIndent();
}

void TBamFile::startProgressReporting(bool indent) const {
	_lastProgressPrinted = 0;
	_timer.start();

	if (indent) logfile().startIndent("Parsing through BAM file ",_filename , ":");
	else logfile().list("Parsing through BAM file ",_filename , ".");
}

void TBamFile::printProgress() const {
	if (atStart() || chrChanged()) {
		logfile().list("Parsing Chromomsome ", curChromosome().name(), ".");
	}
	if (_numAlignmentRead - _lastProgressPrinted >= _progressFrequency) {
		logfile().list("Parsed " + impl::millionReadsRead(_numAlignmentRead) + " million reads (est. " +
					   coretools::str::toStringWithPrecision((100.*_bamReader.Tell())/_fileSize, 2) + "%) in " +
					   _timer.formattedTime());
		_lastProgressPrinted = _numAlignmentRead;
	}
}

void TBamFile::printEndWithSummary(std::string_view outputName, bool indent) const {
	logfile().list("Reached end of BAM file ", _filename, " in " + _timer.formattedTime() + ':');
	logfile().conclude("Parsed a total of " + impl::millionReadsRead(_numAlignmentRead) + " million reads in " + _timer.formattedTime() + '.');

	if (indent) logfile().endIndent();

	printSummary(outputName);
}
void TBamFile::setCurReadGroup(size_t id) {
	if (id != _read.readGroupID) { _read.bamAlignment.EditTag("RG", "Z", _readGroups.getName(id)); }
}

void TBamFile::addCurSamField(const std::string& tag, float value){
       if(_read.bamAlignment.HasTag(tag)){
               _read.bamAlignment.EditTag(tag, "f", value);
       } else {
               _read.bamAlignment.AddTag(tag, "f", value);
       }
}




} //end namespace

