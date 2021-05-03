/*
 * TBamFile.cpp
 *
 *  Created on: May 18, 2020
 *      Author: wegmannd
 */

#include "TBamFile.h"

namespace BAM{

//-----------------------------------------------------
//TBamFile
//-----------------------------------------------------
TBamFile::TBamFile(){
	_open = false;
	_fileSize = 0;
	_numAlignmentRead = 0;
	_numAlignmentsPassedQC = 0;
	_limitNumReads = false;
	_maxNumReadsToRead = 0;

	//current alignment
	_curReadGroupID = 0;
	_chrChanged = false;

	//set filters to default
	_QCFiltersPassed = false;
	_maxReadLength = 65535;
	_keepAll = true; //by default, keep all reads

	//blacklist
	_updateLog = false;
	_bamLog = nullptr;

	//progress reporting
	_logfile = nullptr;
	_progressFrequency = 100000;
	_lastProgressPrinted = 0;
};

void TBamFile::setLimits(TParameters & params, TLog* logfile){
	//number of reads
	if(params.parameterExists("limitReads")){
		_maxNumReadsToRead = params.getParameterInt("limitReads");
		logfile->list("Will limit the analysis to the first " + toString(_maxNumReadsToRead) + " reads in the BAM file.");
		_limitNumReads = true;
	}

	//limit chromosomes?
	_chromosomes.limitAndSetPloidy(params, logfile);

	//limit read groups
	if(params.parameterExists("readGroup")){
		_readGroups.filterReadGroups(params.getParameterString("readGroup"));
		logfile->startIndent("Will limit analysis to the following read groups:");
		_readGroups.printReadgroupsInUse(logfile);
		logfile->endIndent();
		_readGroupFilter.filter("Read group not in use");
	} else {
		_readGroupFilter.keep();
	}
};

void TBamFile::setFilters(TParameters & params, TLog* logfile){
	//max read length
	//is relevant for storage, so only accept values up to 2^16
	//Others will be filtered out.
	int MaxReadLength = params.getParameterIntWithDefault("maxReadLength", 200);
	logfile->list("Expect no read to be longer than " + toString(MaxReadLength) + " bp. (parameter 'maxReadLength')");
	if(MaxReadLength < 1)
		throw "Max read length must be at least 1 bp!";
	if(MaxReadLength > 65536)
		throw "Atlas currently only supports read length up to 65536 bp. Contact us if you have longer reads / fragments";
	_maxReadLength = MaxReadLength;
	_readLengthFilter.filter(1, _maxReadLength, "Read length too large");
	_allowTooLongReads = params.parameterExists("allowTooLongReads");

	//alignment filters
	logfile->startIndent("Will use the following filters on reads:");
	if(params.parameterExists("keepAllReads")){
		_keepAll = true;
		logfile->list("Will keep all reads. (parameter 'keepAllReads', overrules any other QC filter)");
	} else {
		_keepAll = false;
		//duplicates
		if(params.parameterExists("keepDuplicates")){
			_duplicateFilter.keep();
			logfile->list("Duplicate reads: keep. (parameter 'keepDuplicates')");
		} else {
			_duplicateFilter.filter("Duplicate");
			logfile->list("Duplicate reads: filter out. (use 'keepDuplicates' to keep)");
		}

		//soft clips
		if(params.parameterExists("filterSoftClips")){
			_softClippedFilter.filter("Soft clipped");
			logfile->list("Soft clipped reads: filter out. (parameter 'filterSoftClips')");
		} else {
			_softClippedFilter.keep();
			logfile->list("Soft clipped reads: keep. (use 'filterSoftClips' to filter out)");
		}

		//improper pairs
		if(params.parameterExists("keepImproperPairs")){
			_improperPairsFilter.keep();
			logfile->list("Improper pairs: keep. (parameter 'keepImproperPairs')");
		} else {
			_improperPairsFilter.filter("Improper pair");
			logfile->list("Improper pairs: filter out. (use 'keepImproperPairs' to keep)");
		}

		//unmapped reads
		if(params.parameterExists("keepUnmappedReads")){
			_unmappedFilter.keep();
			logfile->list("Unmapped reads: keep. (parameter 'keepUnmappedReads')");
		} else {
			_unmappedFilter.filter("Unmapped");
			logfile->list("Unmapped reads: filter out. (use 'keepUnmappedReads' to keep)");
		}

		//failed QC
		if(params.parameterExists("keepFailedQC")){
			_failedQCFilter.keep();
			logfile->list("Failed QC: keep. (parameter 'keepFailedQC')");
		} else {
			_failedQCFilter.filter("Failed QC");
			logfile->list("Failed QC: filter out. (use 'keepFailedQC' to keep)");
		}

		//secondary reads
		if(params.parameterExists("keepSecondaryReads")){
			_secondaryFilter.keep();
			logfile->list("Secondary reads: keep. (parameter 'keepSecondaryReads')");
		} else {
			_secondaryFilter.filter("Secondary alignment");
			logfile->list("Secondary reads: filter out. (use 'keepSecondaryReads' to keep)");
		}

		//supplementary reads
		if(params.parameterExists("keepSupplementaryReads")){
			_supplementaryFilter.keep();
			logfile->list("Supplementary reads: keep. (parameter 'keepSupplementaryReads')");
		} else {
			_supplementaryFilter.filter("Supplementary alignment");
			logfile->list("Supplementary reads: filter out. (use 'keepSupplementaryReads' to keep)");
		}

		//fragment length
		if(params.parameterExists("keepReadsLongerThanFragment")){
			_longerThanFragmentFilter.keep();
			logfile->list("Reads longer than fragment size: keep. (parameter 'keepReadsLongerThanFragment')");
		} else {
			_longerThanFragmentFilter.filter("Longer than fragment");
			logfile->list("Reads longer than fragment size: filter out. (use 'keepReadsLongerThanFragment' to keep)");
		}

		//strand
		if(params.parameterExists("keepOnlyFwd")){
			_fwdStrandFilter.keep();
			_revStrandFilter.filter("Reverse strand");
			logfile->list("Strand: keep only forward. (parameter 'keepOnlyFwd')");
		}
		else if(params.parameterExists("keepOnlyRev")){
			_fwdStrandFilter.filter("Forward strand");
			_revStrandFilter.keep();
			logfile->list("Strand: keep only reverse. (parameter 'keepOnlyRev')");
		} else {
			_fwdStrandFilter.keep();
			_revStrandFilter.keep();
			logfile->list("Strand: keep forward and reverse. (use 'keepOnlyFwd' or 'keepOnlyRev' to limit)");
		}

		//mate
		if(params.parameterExists("keepOnlyFirst")){
			_firstMateFilter.keep();
			_secondMateFilter.filter("Second mate");
			logfile->list("Mate: keep only first. (parameter 'keepOnlyFirst')");
		}
		else if(params.parameterExists("keepOnlySecond")){
			_firstMateFilter.filter("First mate");
			_secondMateFilter.keep();
			logfile->list("Mate: keep only second. (parameter 'keepOnlySecond')");
		} else {
			_firstMateFilter.keep();
			_secondMateFilter.keep();
			logfile->list("Mate: keep first and second. (use 'keepOnlyFirst' or 'keepOnlySecond' to limit)");
		}

		//blacklist
		if(params.parameterExists("blacklist")){
			std::string blacklistFilename = params.getParameterFilename("blacklist");
			logfile->list("Will filter out reads present in the file '" + blacklistFilename + "'. (parameter 'blacklist')");
			_blacklist.addFromFile(blacklistFilename);
			_blacklistFilter.filter("Was in provided blacklist");
		} else {
			_blacklistFilter.keep();
		}

		//Mapping quality filter
		if(params.parameterExists("minMQ") || params.parameterExists("maxMQ")){
			int MinMQ = params.getParameterIntWithDefault("minMQ", 0);
			int MaxMQ = params.getParameterIntWithDefault("maxMQ", 254);

			if(MinMQ < 0 || MinMQ > 254)
				throw "minMQ '" + toString(MinMQ) + "' is outside the accepted range [0,254]!";
			if(MaxMQ < 0 || MaxMQ > 254)
				throw "maxMQ '" + toString(MaxMQ) + "' is outside the accepted range [0,254]!";
			if(MaxMQ < MinMQ)
				throw "minMQ must be <= maxMQ";

			_mappingQualityFilter.filter(MinMQ, MaxMQ, "Mapping quality outside [" + toString(MinMQ) + ", " + toString(MaxMQ) + "]");
			logfile->list("Mapping quality: restrict to range [" + toString(MinMQ) + ", " + toString(MaxMQ) + "]. (parameters 'minMQ', 'maxMQ')");
		} else {
			_mappingQualityFilter.keep();
			logfile->list("Mapping quality: keep all. (use 'minMQ' and 'maxMQ' to limit)");
		}

		//Fragment length filter
		if(params.parameterExists("minFragmentLength") || params.parameterExists("maxFragmentLength")){
			int MinFragmentLength = params.getParameterIntWithDefault("minFragmentLength", 0);
			int MaxFragmentLength = params.getParameterIntWithDefault("maxFragmentLength", 1000);
			if(MinFragmentLength < 0)
				throw "minFragmentLength '" + toString(MinFragmentLength) + "' must be >0!";
			if(MaxFragmentLength < 0)
				throw "maxFragmentLength '" + toString(MaxFragmentLength) + "' must be >0!";
			if(MinFragmentLength > MaxFragmentLength)
				throw "minFragmentLength must be <= maxFragmentLength!";

			_fragmentLengthFilter.filter(MinFragmentLength, MaxFragmentLength, "Fragment length outside [" + toString(MinFragmentLength) + ", " + toString(MaxFragmentLength) + "]");
			logfile->list("Fragment length: restrict to range [" + toString(MinFragmentLength) + ", " + toString(MaxFragmentLength) + "]. (parameters 'minFragmentLength', 'maxFragmentLength')");
		} else {
			_fragmentLengthFilter.keep();
			logfile->list("Fragment length: keep all. (use 'minFragmentLength', 'maxFragmentLength' to limit)");
		}
	}
	logfile->endIndent();

	//log filtered reads?
	openBamLog(params, logfile);
};

void TBamFile::curFilterOut(){
	_externalFilter.filterOut(_curBamAlignment.Name, _curBamAlignment.IsReverseStrand());
};

void TBamFile::filterOut(const std::string & alignmentName, const bool & isReverseStrand){
	_externalFilter.filterOut(alignmentName, isReverseStrand);
};

void TBamFile::setExternalFilterReason(const std::string reason){
	_externalFilter.setReason(reason);
};

void TBamFile::openBamLog(TParameters & params, TLog* logfile){
	if(params.parameterExists("bamLog") && !_updateLog){
		std::string logFilename = params.getParameterString("bamLog");
		if(logFilename.empty()){
			logFilename = _filename;
			logFilename = extractBeforeLast(logFilename, ".");
			logFilename += ".bamlog.txt.gz";
		}
		logfile->list("Will write all filtered out reads to '" + logFilename + "'.");
		_bamLog = new TBamFileLog(logFilename);
		_updateLog = true;

		//_log to filters
		_duplicateFilter.setLog(_bamLog);
		_softClippedFilter.setLog(_bamLog);
		_improperPairsFilter.setLog(_bamLog);
		_unmappedFilter.setLog(_bamLog);
		_failedQCFilter.setLog(_bamLog);
		_secondaryFilter.setLog(_bamLog);
		_supplementaryFilter.setLog(_bamLog);
		_longerThanFragmentFilter.setLog(_bamLog);
		_readGroupFilter.setLog(_bamLog);
		_fwdStrandFilter.setLog(_bamLog);
		_revStrandFilter.setLog(_bamLog);
		_firstMateFilter.setLog(_bamLog);
		_secondMateFilter.setLog(_bamLog);
		_blacklistFilter.setLog(_bamLog);
		_mappingQualityFilter.setLog(_bamLog);
		_fragmentLengthFilter.setLog(_bamLog);
		_externalFilter.setLog(_bamLog);
	}
};

void TBamFile::writeToBamLog(const std::string & alignmentName, const bool & isReverseStrand, const std::string & reason){
	if(_updateLog){
		_bamLog->write(alignmentName, isReverseStrand, reason);
	}
};

void TBamFile::_fillSamHeader(TSamHeader & SamHeader){
	//Note: chromosomes and read groups are in separate objects
	SamHeader.set(_bamHeader.Version, _bamHeader.SortOrder, _bamHeader.GroupOrder, "none");

	//add programs
	for(auto it = _bamHeader.Programs.Begin(); it != _bamHeader.Programs.End(); ++it){
		SamHeader.addProgram(it->ID, it->Name, it->CommandLine, "", it->Version);
	}

	//add links among programs
	for(auto it = _bamHeader.Programs.Begin(); it != _bamHeader.Programs.End(); ++it){
		if(it->HasPreviousProgramID()){
			SamHeader.addPreviousProgramInChain(it->ID, it->PreviousProgramID);
		}
	}

	//add comments
	for(auto& c : _bamHeader.Comments){
		SamHeader.addComment(c);
	}
};

void TBamFile::_fillChromosomes(TChromosomes & Chromosomes){
	if(_bamHeader.Sequences.Size() < 1){
		throw "No chromosomes present in BAM header!";
	}

	//make sure object is empty
	Chromosomes.clear();

	//copy from BamHeader
	for(BamTools::SamSequenceIterator chrIt=_bamHeader.Sequences.Begin(); chrIt!=_bamHeader.Sequences.End(); ++chrIt){
		Chromosomes.appendChromosome(chrIt->Name, convertString<uint64_t>(chrIt->Length));
	}
};

void TBamFile::_fillReadGroups(TReadGroups & ReadGroups){
	//make sure they are empty
	ReadGroups.clear();

	//now add one by one
	for(auto it = _bamHeader.ReadGroups.Begin(); it != _bamHeader.ReadGroups.End(); ++it){
		//add read group
		const TReadGroup& rg = ReadGroups.add(it->ID);

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
void TBamFile::open(const std::string Filename, const bool IndexNotRequired, TLog* Logfile){
	_logfile = Logfile;
	_filename = Filename;

	//open BAM file
	_logfile->listFlushDots("Opening BAM file '" + _filename);
	if (!_bamReader.Open(_filename))
		throw "Failed to open BAM file '" + Filename + "'!";

	//load index file
	if(!_bamReader.LocateIndex() && !IndexNotRequired)
		throw "No index file found for BAM file '" + Filename + "'!";
	_open = true;

	//initialize bam stuff
	_bamHeader = _bamReader.GetHeader();

	_fillSamHeader(_samHeader);

	//initialize read groups
	_fillReadGroups(_readGroups);

	//initialize chromosomes and set cur chromosome to end
	_fillChromosomes(_chromosomes);
	_curChromosome = _chromosomes.end();

	//get file size
	_bamReader.Jump(_chromosomes.size() - 1, 0);
	BamTools::BamAlignment bamAlignment;
	_bamReader.GetNextAlignment(bamAlignment);
	_fileSize = _bamReader.tell();
	_bamReader.Rewind();

	_logfile->done();
};

void TBamFile::close(){
	if(_open){
		_bamReader.Close();
		_open = false;
	}
};

void TBamFile::_applyFilters(){
	//check read length
	//read length is special as it affects our storage
	if(!_readLengthFilter.pass(_curBamAlignment.AlignedBases.size(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())){
		if(!_allowTooLongReads){
			throw "Alignment '" +  _curBamAlignment.Name + "' is longer than the max read length " + toString(_maxReadLength) + "!\nPlease change max read length to parse this data, or add allowTooLongReads to ignore this error and filter out such reads.";
		} else {
			_QCFiltersPassed = false;
		}
	} else if(_keepAll){
		//keep all?
		_QCFiltersPassed = true;
	} else {
		//apply regular filters
		_QCFiltersPassed =  _duplicateFilter.pass(_curBamAlignment.IsDuplicate(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _softClippedFilter.pass(_curCigar.lengthSoftClipped() > 0, _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _improperPairsFilter.pass(_curBamAlignment.IsPaired() && !_curBamAlignment.IsProperPair(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _unmappedFilter.pass(!_curBamAlignment.IsMapped(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _failedQCFilter.pass(_curBamAlignment.IsFailedQC(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _secondaryFilter.pass(!_curBamAlignment.IsPrimaryAlignment(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _supplementaryFilter.pass(_curBamAlignment.IsSupplementary(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _readGroupFilter.pass(!_readGroups.readGroupInUse(_curReadGroupID), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _fwdStrandFilter.pass(!_curBamAlignment.IsReverseStrand(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _revStrandFilter.pass(_curBamAlignment.IsReverseStrand(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _firstMateFilter.pass(_curBamAlignment.IsFirstMate(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _secondMateFilter.pass(_curBamAlignment.IsSecondMate(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _mappingQualityFilter.pass(_curBamAlignment.MapQuality, _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _blacklistFilter.pass(_blacklist.isInBlacklist(_curBamAlignment.Name), _curBamAlignment.Name, _curBamAlignment.IsSecondMate());

		//fragment length
		if(_QCFiltersPassed){
			_QCFiltersPassed = _fragmentLengthFilter.pass(curFragmentLength(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
                               && _longerThanFragmentFilter.pass(_curBamAlignment.IsProperPair() && _curBamAlignment.InsertSize < _curCigar.lengthAligned(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate());
		}
	}

	//update counter
	if(_QCFiltersPassed){
		++_numAlignmentsPassedQC;
	}
};

bool TBamFile::readNextAlignment(){
	//check if we limit reads
	if(_limitNumReads && _numAlignmentRead >=_maxNumReadsToRead){
		return false;
	}

	//store previous position
	_previousAlignmentPosition = _curAlignmentPosition;

	//get next alignment
	if(!_bamReader.GetNextAlignment(_curBamAlignment)){
		return false;
	}

	//check if chromosome changed
	if(_curChromosome == _chromosomes.end() || _curBamAlignment.RefID != _curChromosome->refID()){
		//advance chromosome
		if(_curChromosome == _chromosomes.end()){
			_curChromosome = _chromosomes.begin();
		}

		while(_curBamAlignment.RefID != _curChromosome->refID()){
			++_curChromosome;

			if(_curChromosome == _chromosomes.end()){
				//is chromosome not in header?
				if(!_chromosomes.exists(_curBamAlignment.RefID)){
					throw "Chromosome with refID " + toString(_curBamAlignment.RefID) + " is missing from BAM header!";
				} else {
					throw "BAM file not sorted!";
				}
			}
		}

		//if not in use: jump to next in use
		if(!_curChromosome->inUse){
			while(!_curChromosome->inUse){
				++_curChromosome;

				if(_curChromosome == _chromosomes.end()){
					return false;
				}
			}

			//jump reader and read first alignment
			jump(_curChromosome->chrStart);
			if(!_bamReader.GetNextAlignment(_curBamAlignment)){
				return false;
			}
		}
		_chrChanged = true;
	} else {
		_chrChanged = false;
	}

	//get current position and update counter
	_curAlignmentPosition.move(_curBamAlignment.RefID, _curBamAlignment.Position);
	++_numAlignmentRead;


	//check if BAM file is sorted
	if(_curAlignmentPosition < _previousAlignmentPosition){
		throw "BAM file must be sorted by position! Alignment '" + _curBamAlignment.Name + "' is at position " + toString(_curBamAlignment.Position) + ", which is before the position of the previous alignment (" + toString(_previousAlignmentPosition.position()) + ")";
	}

	//store current read group ID
	std::string readGroup;
	_curBamAlignment.GetTag("RG", readGroup);
	_curReadGroupID = _readGroups.getId(readGroup);

	//parse CIGAR
	_curCigar.clear();
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

bool TBamFile::jump(const TGenomePosition Position){
	_previousAlignmentPosition.clear();
	return _bamReader.Jump(Position.refID(), Position.position());
};

void TBamFile::rewind(){
	_bamReader.Rewind();
	_numAlignmentRead = 0;
	_numAlignmentsPassedQC = 0;
	_previousAlignmentPosition.clear();
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
		throw "Failed to open BAM file '" + filename + "'!";
};

void TBamFile::writeCurAlignment(TOutputBamFile & out){
	out._writeAlignment(_curBamAlignment);
};

//--------------------------------------------------------
// Getters and setters of cur alignment
//--------------------------------------------------------
uint16_t TBamFile::curFragmentLength() const{
	if(_curBamAlignment.IsProperPair()){
		return abs(_curBamAlignment.InsertSize) + _curCigar.lengthInserted() - _curCigar.lengthDeleted();
	} else {
		return _curCigar.lengthSequenced();
	}
};

uint16_t TBamFile::curUsableAlignedLength(TQualityFilter & qualFilter) const{
	uint16_t counter = 0;
	for(size_t d=0; d<_curBamAlignment.AlignedBases.length(); ++d){
		if(_curBamAlignment.AlignedBases.at(d) != N && qualFilter.pass(_curBamAlignment.AlignedQualities.at(d))){
			++counter;
		}
	}
	return counter;
};

std::string TBamFile::curQuerySequence(const uint16_t start, const uint16_t length) const{
	return _curBamAlignment.QueryBases.substr(start, length);
};

void TBamFile::curSetNewReadGroup(const uint16_t id){
	if(id != _curReadGroupID){
		_curBamAlignment.EditTag("RG", "Z", _readGroups.getName(id));
	}
};

void TBamFile::curAddSamField(const std::string tag, const std::string value){
	if(_curBamAlignment.HasTag(tag)){
		_curBamAlignment.EditTag(tag, "Z", value);
	} else {
		_curBamAlignment.AddTag(tag, "Z", value);
	}
};

void TBamFile::curAddSamField(const std::string tag, const float value){
	if(_curBamAlignment.HasTag(tag)){
		_curBamAlignment.EditTag(tag, "f", value);
	} else {
		_curBamAlignment.AddTag(tag, "f", value);
	}
};

//-----------------------------------------------------
// Reporting
//-----------------------------------------------------
void TBamFile::printSummaryNoEndIndent(){
	_logfile->startIndent("Summary of parsed reads from BAM file '" + _filename + "':");
	_logfile->list("Total number of reads read: " + toString(_numAlignmentRead));
	_logfile->list("Reads that passed filters: " + toString(_numAlignmentsPassedQC) + " (" + toPercentString(_numAlignmentsPassedQC, _numAlignmentRead, 3) + "%)");
	uint64_t numFiltered = _numAlignmentRead - _numAlignmentsPassedQC;
	_logfile->list("Reads that were filtered out: " + toString(numFiltered) + " (" + toPercentString(numFiltered, _numAlignmentRead, 3) + "%)");

	if(numFiltered > 0){
		_logfile->addIndent();

		_duplicateFilter.summary(_logfile, numFiltered);
		_softClippedFilter.summary(_logfile, numFiltered);
		_improperPairsFilter.summary(_logfile, numFiltered);
		_unmappedFilter.summary(_logfile, numFiltered);
		_failedQCFilter.summary(_logfile, numFiltered);
		_secondaryFilter.summary(_logfile, numFiltered);
		_supplementaryFilter.summary(_logfile, numFiltered);
		_longerThanFragmentFilter.summary(_logfile, numFiltered);
		_readGroupFilter.summary(_logfile, numFiltered);
		_fwdStrandFilter.summary(_logfile, numFiltered);
		_revStrandFilter.summary(_logfile, numFiltered);
		_firstMateFilter.summary(_logfile, numFiltered);
		_secondMateFilter.summary(_logfile, numFiltered);
		_blacklistFilter.summary(_logfile, numFiltered);
		_mappingQualityFilter.summary(_logfile, numFiltered);
		_fragmentLengthFilter.summary(_logfile, numFiltered);
		_externalFilter.summary(_logfile, numFiltered);
		_logfile->endIndent();
	}
};

void TBamFile::printSummary(){
	printSummaryNoEndIndent();
	_logfile->endIndent();
	_logfile->endIndent();
};

void TBamFile::startProgressReporting(uint32_t Frequency){
	if(!_open){
		throw "Can not start progress reporting of BAM file: BAM file not open!";
	}

	_progressFrequency = Frequency;
	_lastProgressPrinted = 0;
	_timer.start();

	_logfile->startIndent("Parsing through BAM file:");
};

void TBamFile::printProgress(){
	if(_numAlignmentRead - _lastProgressPrinted >= _progressFrequency){
		_logfile->list("Parsed " + _millionReadsRead() + " million reads (est. " + to_string_with_precision(positionInFile() * 100, 2) + "%) in " + toString(_timer.minutes()) + " min.");
		_lastProgressPrinted = _numAlignmentRead;
	}
};

void TBamFile::printEndWithSummary(){
	printEndNoEndIndent();
	_logfile->endIndent();
	printSummary();
};

void TBamFile::printEndNoEndIndent(){
	_logfile->list("Reached end of BAM file in " + toString(_timer.minutes()) + " min.");
	_logfile->conclude("Parsed a total of " + _millionReadsRead() + " million reads in " + toString(_timer.minutes()) + " min.");
};

//-----------------------------------------------------
//TOutputBamFile
//----------------------------------------------------
TOutputBamFile::TOutputBamFile(){
	_openForWriting = false;
	_readGroups = nullptr;
	_genoMap = nullptr;
	_qualityMap = nullptr;
};

TOutputBamFile::TOutputBamFile(const std::string Filename, const TBamFile & Original, GenotypeLikelihoods::TGenotypeMap* GenoMap, TQualityMap* QualityMap){
	_openForWriting = false;
	open(Filename, Original.samHeader(), Original.chromosomes(), Original.readGroups(), GenoMap, QualityMap);
};

TOutputBamFile::~TOutputBamFile(){
	closeNoIndex();
};

void TOutputBamFile::open(const std::string Filename, const TSamHeader & Header, const TChromosomes & Chromosomes, const TReadGroups & ReadGroups, GenotypeLikelihoods::TGenotypeMap* GenoMap, TQualityMap* QualityMap){
	closeNoIndex();

	_outputFilename = Filename;
	_readGroups = &ReadGroups;
	_genoMap = GenoMap;
	_qualityMap = QualityMap;

	//construct new header /without chromosomes
	std::string header = Header.compileSamHeader(ReadGroups, Chromosomes);

	//fill bamtools chromosomes
	BamTools::RefVector ref;
	for(auto it = Chromosomes.cbegin(); it != Chromosomes.cend(); ++it){
		ref.emplace_back(it->name, it->length);
	}

	//open file for writing
	if(!_bamWriter.Open(_outputFilename, header, ref)){
		throw "Failed to open BAM file '" + _outputFilename + "'!";
	}

	_openForWriting = true;
};

void TOutputBamFile::open(const std::string Filename, const TBamFile & Original, GenotypeLikelihoods::TGenotypeMap* GenoMap, TQualityMap* QualityMap){
	open(Filename, Original.samHeader(), Original.chromosomes(), Original.readGroups(), GenoMap, QualityMap);
};

void TOutputBamFile::close(TLog* logfile){
	if(_openForWriting){
		_bamWriter.Close();

		logfile->listFlush("Creating index of BAM file '" + _outputFilename + "' ...");
		BamTools::BamReader reader;
		if(!reader.Open(_outputFilename))
			throw "Failed to open BAM file '" + _outputFilename + "' for indexing!";

		// create index for BAM file
		reader.CreateIndex(BamTools::BamIndex::STANDARD);

		//close BAM file
		reader.Close();
		logfile->done();

		_openForWriting = false;
	}
};

void TOutputBamFile::close(){
	if(_openForWriting){
		_bamWriter.Close();

		// create index of BAM file
		BamTools::BamReader reader;
		if(!reader.Open(_outputFilename))
            throw "Failed to open BAM file '" + _outputFilename + "' for indexing!";
		reader.CreateIndex(BamTools::BamIndex::STANDARD);

        //close BAM file
		reader.Close();
		_openForWriting = false;
	}
};

void TOutputBamFile::closeNoIndex(){
	if(_openForWriting){
		_bamWriter.Close();
		_openForWriting = false;
	}
};

void TOutputBamFile::_writeAlignment(BamTools::BamAlignment & alignment){
	if(!_openForWriting){
		throw "BAM writer is not open!";
	}

	//adjust qualities for printing
	_qualityMap->adjustQualitiesForWriting(alignment.Qualities);

	// write alignment
	if(!_bamWriter.SaveAlignment(alignment))
		throw "Read '" + alignment.Name + "' could not be written!";
};

void TOutputBamFile::_writeAlignment(const TAlignment & alignment){
	if(!_openForWriting){
		throw "BAM writer is not open!";
	}

	//create bamAlignment and then write
	BamTools::BamAlignment _tmpBamAlignment;

	_tmpBamAlignment.Name = alignment.name();
	_tmpBamAlignment.AlignmentFlag = alignment.flags();
	_tmpBamAlignment.RefID = alignment.refID();
	_tmpBamAlignment.Position = alignment.position();
	_tmpBamAlignment.InsertSize = alignment.insertSize();
	_tmpBamAlignment.MapQuality = alignment.mappingQuality();

	if(alignment.isPaired()){
		_tmpBamAlignment.MateRefID = alignment.mateRefID();
		_tmpBamAlignment.MatePosition = alignment.matePosition();
	}

	//CIGAR
	for(const auto& it : alignment.cigar()){
		_tmpBamAlignment.CigarData.emplace_back(it.type, it.length);
	}

	//add sequences and qualities
	_tmpBamAlignment.QueryBases = alignment.sequence(*_genoMap, *_qualityMap);
	_tmpBamAlignment.Qualities = alignment.qualities(*_genoMap, *_qualityMap);

	//add read group information
	_tmpBamAlignment.AddTag("RG", "Z", _readGroups->getName(alignment.readGroupId()));

	//and now write
	if(!_bamWriter.SaveAlignment(_tmpBamAlignment))
		throw "Read '" + _tmpBamAlignment.Name + "' could not be written!";
};

void TOutputBamFile::writeAlignment(const TAlignment & alignment){
	auto it = _futureAlignments.begin();
	if(it != _futureAlignments.end()){
		//write alignments BEFORE alignment to write next
		while(it != _futureAlignments.end() && *it < alignment){
			_writeAlignment(*it);
		}

		//remove written alignments
		_futureAlignments.erase(_futureAlignments.begin(), it);
	}

	//write next alignment
	_writeAlignment(alignment);
};

}; //end namespace

