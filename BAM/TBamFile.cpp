/*
 * TBamFile.cpp
 *
 *  Created on: May 18, 2020
 *      Author: wegmannd
 */

#include "TBamFile.h"

#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <exception>
#include <filesystem>
#include <string>

#include "api/BamIndex.h"
#include "api/SamProgram.h"
#include "api/SamProgramChain.h"
#include "api/SamReadGroup.h"
#include "api/SamReadGroupDictionary.h"
#include "api/SamSequence.h"
#include "api/SamSequenceDictionary.h"

#include "coretools/Main/TLog.h"
#include "coretools/Math/TNumericRange.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/globalConstants.h"
#include "coretools/Types/strongTypes.h"

#include "TSamFlags.h"
#include "TBamFilters.h"
#include "TOutputBamFile.h"

namespace BAM{
using coretools::instances::parameters;
using coretools::instances::logfile;

void TBamFile::_setLimits(){
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
	}
};

void TBamFile::curFilterOut(){
	_filters.filterOut(FilterType::External, _curBamAlignment.Name, _curBamAlignment.IsReverseStrand(), _curReadGroupID, refID());
};

void TBamFile::filterOut(const TAlignment & Alignment){
	_filters.filterOut(FilterType::External, Alignment.name(), Alignment.isReverseStrand(), Alignment.readGroupId(), Alignment.refID());

};

void TBamFile::setExternalFilterReason(std::string_view reason){
	_filters.enable(FilterType::External, reason);
};

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
};

void TBamFile::_fillChromosomes(){
	if(_bamHeader.Sequences.Size() < 1){
		UERROR("No chromosomes present in BAM header!");
	}

	//make sure object is empty
	_chromosomes.clear();

	//copy from BamHeader
	for(BamTools::SamSequenceIterator chrIt=_bamHeader.Sequences.Begin(); chrIt!=_bamHeader.Sequences.End(); ++chrIt){
		_chromosomes.appendChromosome(chrIt->Name, coretools::str::fromString<uint64_t>(chrIt->Length));
	}
};

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
};

//--------------------------------------------------------
// Functions for reading
//--------------------------------------------------------
TBamFile::TBamFile(std::string_view Filename){
	_filename = Filename;

	//open BAM file
	logfile().list("Opening BAM file '", _filename, "'.");
	if (!_bamReader.Open(_filename))
		UERROR("Failed to open BAM file '", Filename, "'!");

	//load or create index file
	const std::string fnIndex1 = std::string(Filename).append(".bai");
	Filename.remove_suffix(4);
	const std::string fnIndex2 = std::string(Filename).append(".bai");
	if (std::filesystem::exists(fnIndex1)) {
		logfile().list("Opening BAM index file '", fnIndex1, "'.");
		if(!_bamReader.OpenIndex(fnIndex1))
			UERROR("Failed to open BAM index file '", fnIndex1, "'!");
	}
	else if (std::filesystem::exists(fnIndex2)) {
		logfile().list("Opening BAM index file '", fnIndex2, "'.");
		if (!_bamReader.OpenIndex(fnIndex2))
			UERROR("Failed to open BAM index file '", fnIndex2, "'!");
	} else {
		logfile().list("Creating BAM index file '", fnIndex1, "'.");
		if (!_bamReader.CreateIndex())
			UERROR("Failed to create BAM index file '", fnIndex1, "'!");
	}

	//initialize bam stuff
	_bamHeader = _bamReader.GetHeader();

	_fillSamHeader();

	//initialize read groups
	_fillReadGroups();
	_numNotAligned.resize(_readGroups.size());

	//initialize chromosomes and set cur chromosome to end
	_fillChromosomes();
	_curChromosome = _chromosomes.end();

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

	_filters.resize(_readGroups.size(), _chromosomes.size(), _filename);
	_setLimits();
};

void TBamFile::setFilters(const TBamFilters& Filters) {
	_filters.clone(Filters);
}

bool TBamFile::_readNextAlignmentFromFile(){
	if(!_bamReader.GetNextAlignment(_curBamAlignment)){
		return false;
	}
	++_numAlignmentRead;

	//store current read group ID
	std::string readGroup;
	_curBamAlignment.GetTag("RG", readGroup);
	_curReadGroupID = _readGroups.getId(readGroup);

	return true;
}

void TBamFile::_applyFilters() {
	// MappedLength filter is always set
	if (!_filters.pass(FilterType::MappedLength, _curCigar.lengthMapped(), _curBamAlignment.Name,
					   _curBamAlignment.IsSecondMate(), _curReadGroupID, refID())) {
		_QCFiltersPassed = false;
	} else if (!_filters.enabled()) {
		_QCFiltersPassed = true;
	} else {
	// apply regular filters
	_QCFiltersPassed =
		_filters.pass(FilterType::Duplicate, !_curBamAlignment.IsDuplicate(), _curBamAlignment.Name,
		              _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::SoftClippedRation,
		              static_cast<double>(_curCigar.lengthSoftClipped()) / _curCigar.lengthRead() <=
					  _filters.softClipRation(),
		              _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::ImproperPairs, !_curBamAlignment.IsPaired() || _curBamAlignment.IsProperPair(),
		              _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::Unmapped, _curBamAlignment.IsMapped(), _curBamAlignment.Name,
		              _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::FailedQC, !_curBamAlignment.IsFailedQC(), _curBamAlignment.Name,
		              _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::Secondary, _curBamAlignment.IsPrimaryAlignment(), _curBamAlignment.Name,
		              _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::Supplementary, !_curBamAlignment.IsSupplementary(), _curBamAlignment.Name,
		              _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::ReadGroup, _readGroups.readGroupInUse(_curReadGroupID), _curBamAlignment.Name,
		              _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::FwdStrand, _curBamAlignment.IsReverseStrand(), _curBamAlignment.Name,
		              _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::RevStrand, !_curBamAlignment.IsReverseStrand(), _curBamAlignment.Name,
		              _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::FirstMate, _curBamAlignment.IsFirstMate(), _curBamAlignment.Name,
		              _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::SecondMate, _curBamAlignment.IsSecondMate(), _curBamAlignment.Name,
		              _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::MappingQuality, (size_t)_curBamAlignment.MapQuality, _curBamAlignment.Name,
		              _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::Blacklist, !_filters.blacklist().isInBlacklist(_curBamAlignment.Name), _curBamAlignment.Name,
		              _curBamAlignment.IsSecondMate(), _curReadGroupID, refID()) &&
		_filters.pass(FilterType::ReadLength, _curCigar.lengthRead(), _curBamAlignment.Name,
		              _curBamAlignment.IsSecondMate(), _curReadGroupID, refID());

	// fragment length
	if (_QCFiltersPassed) {
			_QCFiltersPassed = _filters.pass(FilterType::FragmentLength, curFragmentLength(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID, refID())
				&& _filters.pass(FilterType::LongerThanFragment, !_curBamAlignment.IsProperPair() || abs(_curBamAlignment.InsertSize) >= static_cast<int32_t>(_curCigar.lengthAligned()), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID, refID());
		}
	}

	//update counter
	if(_QCFiltersPassed){
		++_numAlignmentsPassedQC;
	}
};

bool TBamFile::readNextAlignment(){
	//check if we limit reads
	if(_numAlignmentRead >=_maxNumReadsToRead){
		return false;
	}

	//store previous position
	_previousAlignmentPosition = _curAlignmentPosition;


	//check if it has no read group
	bool pass = true;
	do {
		// get next alignment
		if (!_readNextAlignmentFromFile()) { return false; }

		if (_curReadGroupID == TReadGroups::noReadGroupId) {
			++_numNoReadGroup;
			pass = false;
		}

		// check if it is unaligned (refID < 0), in which case we read until the first aligned read
		if (_curBamAlignment.RefID < 0) {
			++_numNotAligned[_curReadGroupID];
			pass = false;
		}
	} while (!pass);

	//check if chromosome changed
	if(_curChromosome == _chromosomes.end() || _curBamAlignment.RefID != static_cast<int>(_curChromosome->refID())){
		//advance chromosome
		if(_curChromosome == _chromosomes.end()){
			_curChromosome = _chromosomes.begin();
		}

		while(_curBamAlignment.RefID != static_cast<int>(_curChromosome->refID())){
			++_curChromosome;

			if(_curChromosome == _chromosomes.end()){
				//is chromosome not in header?
				if(!_chromosomes.exists(_curBamAlignment.RefID)){
					UERROR("Chromosome with refID ", _curBamAlignment.RefID, " is missing from BAM header!");
				} else {
					UERROR("BAM file not sorted!");
				}
			}
		}

		//if not in use: jump to next in use
		if(!_curChromosome->inUse()){
			while(!_curChromosome->inUse()){
				++_curChromosome;

				if(_curChromosome == _chromosomes.end()){
					return false;
				}
			}

			//jump reader and read first alignment
			jump(_curChromosome->start());
			if(!_bamReader.GetNextAlignment(_curBamAlignment)){
				return false;
			}
		}
		_chrChanged = true;
	} else {
		_chrChanged = false;
	}

	//get current position, clear CIGAR and update counter
	_curAlignmentPosition.move(_curBamAlignment.RefID, _curBamAlignment.Position);
	_curCigar.clear(); //needs to be cleared here to be empty in case of alignments that are unaligned

	//check if BAM file is sorted
	if(_curAlignmentPosition < _previousAlignmentPosition){
		UERROR("BAM file must be sorted by position! Alignment '", _curBamAlignment.Name, "' is at position ", _curBamAlignment.Position, ", which is before the position of the previous alignment (", _previousAlignmentPosition.position(), ")");
	}

	//update per read group counter
	if(_curReadGroupID != TReadGroups::noReadGroupId){
		_numAlignmentReadPerReadGroupPerChromosome.add(_curReadGroupID, _curChromosome->refID());
	}

	//parse CIGAR
	for(auto& it : _curBamAlignment.CigarData){
		_curCigar.add(it.Type, it.Length);
	}
	//apply filters
	_applyFilters();

	return true;
};

bool TBamFile::readNextAlignmentThatPassesFilters(){
	_QCFiltersPassed = false;
	while(!_QCFiltersPassed){
		if(!readNextAlignment()){
			return false;
		}
	}
	return true;
};

void TBamFile::fill(TAlignment & alignment) const{
	alignment.fill(_curBamAlignment.Name,
				   TSamFlags(_curBamAlignment.AlignmentFlag),
				   _curBamAlignment.RefID,
				   _curBamAlignment.Position,
				   _curBamAlignment.MapQuality,
				   _curCigar,
				   _curBamAlignment.MateRefID,
				   _curBamAlignment.MatePosition,
				   _curBamAlignment.InsertSize,
				   _curBamAlignment.QueryBases,
				   _curBamAlignment.Qualities,
				   _curReadGroupID);
};

bool TBamFile::readNextAlignment(TAlignment & alignment){
	if(!readNextAlignment()){
		return false;
	}

	fill(alignment);
	return true;
};

bool TBamFile::readNextAlignmentThatPassesFilters(TAlignment & alignment){
	if(!readNextAlignmentThatPassesFilters()){
		return false;
	}
	fill(alignment);
	return true;
};

bool TBamFile::jump(const genometools::TGenomePosition Position){
	_previousAlignmentPosition.clear();
	_curAlignmentPosition.clear();
	return _bamReader.Jump(Position.refID(), Position.position());
};

//--------------------------------------------------------
// Functions for writing
//--------------------------------------------------------
void TBamFile::_openForWriting(BamTools::BamWriter & bamWriter, const std::string filename){
	//construct new header
	BamTools::SamHeader newHeader(_bamHeader);

	//make sure read groups are OK: copy from readGroup object
	newHeader.ReadGroups.Clear();
	for(auto it = _readGroups.begin(); it!=_readGroups.end(); ++it){
		if(it->writeToHeader){
			BamTools::SamReadGroup newRg(it->name_ID);

			//copy rest
			newRg.Description = it->description_DS;
			newRg.FlowOrder = it->flowOrder_FO;
			newRg.KeySequence = it->keySequence_KS;
			newRg.Library = it->library_LB;
			newRg.PlatformUnit = it->platformUnit_PU;
			newRg.PredictedInsertSize = it->predictedInsertSize_PI;
			newRg.ProductionDate = it->productionDate_DT;
			newRg.Program = it->program_PG;
			newRg.Sample = it->sample_SM;
			newRg.SequencingCenter = it->sequencingCenter_CN;
			newRg.SequencingTechnology = it->sequencingTechnology_PL;

			//add to header
			newHeader.ReadGroups.Add(newRg);
		}
	}

	//extract references
	BamTools::RefVector references = _bamReader.GetReferenceData();

	//open file for writing
	if(!bamWriter.Open(filename, _bamHeader, references))
		UERROR("Failed to open BAM file '", filename, "'!");
};

void TBamFile::writeCurAlignment(TOutputBamFile & out){
	out.writeAlignment(_curBamAlignment);
};

//--------------------------------------------------------
// Getters and setters of cur alignment
//--------------------------------------------------------
size_t TBamFile::curFragmentLength() const{
	if(_curBamAlignment.IsProperPair()){
		const size_t inserted = abs(_curBamAlignment.InsertSize) + _curCigar.lengthInserted();
		const size_t deleted  = _curCigar.lengthDeleted();
		if (inserted < deleted) return 0;
		return inserted - deleted;
	} else {
		return _curCigar.lengthSequenced();
	}
};

std::string TBamFile::curQuerySequence(size_t start, size_t length) const{
	return _curBamAlignment.QueryBases.substr(start, length);
};

void TBamFile::curSetNewReadGroup(size_t id){
	if(id != _curReadGroupID){
		_curBamAlignment.EditTag("RG", "Z", _readGroups.getName(id));
	}
};

void TBamFile::curAddSamField(const std::string& tag, float value){
	if(_curBamAlignment.HasTag(tag)){
		_curBamAlignment.EditTag(tag, "f", value);
	} else {
		_curBamAlignment.AddTag(tag, "f", value);
	}
};

//-----------------------------------------------------
// Reporting
//-----------------------------------------------------
void TBamFile::_writeFilteringStats(std::string_view outputName) const {
	std::string filename = std::string(outputName).append("_filterSummary.txt");
	coretools::instances::logfile().listFlush("Writing general filter counts to '" + filename + "' ...");

	//creating header
	std::vector<std::string> header;
	header.push_back("readGroup");
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
	logfile().list("Reads that passed filters: ", _numAlignmentsPassedQC, " (", coretools::str::toPercentString(_numAlignmentsPassedQC, _numAlignmentRead, 3), "%)");
	const auto numFiltered = _numAlignmentRead - _numAlignmentsPassedQC;
	logfile().list("Reads that were filtered out: ", numFiltered, " (" + coretools::str::toPercentString(numFiltered, _numAlignmentRead, 3), "%)");

	//write counts of filtered reads for each read group to _filterSummary.txt file
	_writeFilteringStats(outputName);

	//print counts of filtered reads for each read group to terminal, doesn't list filters if no reads were removed
	for (size_t rg = 0; rg < _readGroups.size(); rg++){
		//logfile().newLine();
		logfile().list("Number of reads filtered from read group: '" + coretools::str::toString(_readGroups.getName(rg)) + "'");
		logfile().addIndent();
		if (rg < _numNotAligned.size() && _numNotAligned[rg] > 0) logfile().list("Not aligned: ", _numNotAligned[rg]);
		_filters.summary(numFiltered, rg);
		logfile().endIndent();
	}

	logfile().endIndent();
	logfile().endIndent();
};

void TBamFile::startProgressReporting(size_t Frequency) const {
	_progressFrequency = Frequency;
	_lastProgressPrinted = 0;
	_timer.start();

	logfile().startIndent("Parsing through BAM file:");
};

void TBamFile::printProgress() const {
	if (_numAlignmentRead - _lastProgressPrinted >= _progressFrequency) {
		logfile().list("Parsed " + _millionReadsRead() + " million reads (est. " +
					   coretools::str::toStringWithPrecision((100.*_bamReader.Tell())/_fileSize, 2) + "%) in " +
					   _timer.formattedTime());
		_lastProgressPrinted = _numAlignmentRead;
	}
};

void TBamFile::printEndWithSummary(std::string_view outputName) const {
	logfile().list("Reached end of BAM file in " + _timer.formattedTime() + ':');
	logfile().conclude("Parsed a total of " + _millionReadsRead() + " million reads in " + _timer.formattedTime() + '.');
	logfile().endIndent();
	printSummary(outputName);
};



}; //end namespace

