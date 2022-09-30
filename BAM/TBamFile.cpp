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
#include "TLog.h"
#include "TNumericRange.h"
#include "TParameters.h"
#include "TSamFlags.h"
#include "api/BamIndex.h"
#include "api/SamProgram.h"
#include "api/SamProgramChain.h"
#include "api/SamReadGroup.h"
#include "api/SamReadGroupDictionary.h"
#include "api/SamSequence.h"
#include "api/SamSequenceDictionary.h"
#include "globalConstants.h"
#include "strongTypes.h"

namespace BAM{
using coretools::TParameters;
using coretools::TLog;
using coretools::TNumericRange;
using genometools::BaseQuality;

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
	_keepAll = true; //by default, keep all reads
	_allowTooLongReads = false;

	//progress reporting
	_logfile = nullptr;
	_progressFrequency = 100000;
	_lastProgressPrinted = 0;
};

void TBamFile::setLimits(TParameters & params, TLog* logfile){
	//number of reads
	if(params.parameterExists("limitReads")){
		_maxNumReadsToRead = params.getParameter<uint64_t>("limitReads");
		logfile->list("Will limit the analysis to the first ", _maxNumReadsToRead, " reads in the BAM file.");
		_limitNumReads = true;
	}

	//limit chromosomes?
	_chromosomes.limitAndSetPloidy(params, logfile);

	//limit read groups
	if(params.parameterExists("readGroup")){
		_readGroups.filterReadGroups(params.getParameter<std::string>("readGroup"));
		logfile->startIndent("Will limit analysis to the following read groups:");
		_readGroups.printReadgroupsInUse(logfile);
		logfile->endIndent();
		_readGroupFilter.filter("Read group not in use");
	} else {
		_readGroupFilter.keep();
	}
};

void TBamFile::setFilters(TParameters & params, TLog* logfile){
	//alignment filters
	logfile->startIndent("Will use the following filters on reads:");

	//mapping length
	//--------------
	//is relevant for storage
	//print error if reads are longer and filter is default
	TNumericRange<uint32_t> mappingLengthRange;
	if(params.parameterExists("filterMappingLength")){
		params.fillParameter("filterMappingLength", mappingLengthRange);
		_allowTooLongReads = true;
	} else {
		//set default
		mappingLengthRange.set(0, true, 200, true);
		_allowTooLongReads = params.parameterExists("allowTooLongReads");
	}
	_mappedLengthFilter.filter(mappingLengthRange, "Mapped length outside " + mappingLengthRange.rangeString());
	logfile->list("Mapped length: restrict to range " + _mappedLengthFilter.rangeString() + ". (parameter 'filterMappingLength')");
	if(mappingLengthRange.max() > 100000){
		logfile->warning("The chosen mapping length filter allows for reads to span >100kb of the reference genome. This may affect performance in case of paired-end reads.");
	}

	//keep all otherwise?
	//-------------------
	if(params.parameterExists("keepAllReads")){
		_keepAll = true;
		logfile->list("Will keep all reads. (parameter 'keepAllReads', overrules any other QC filter except filterMappingLength)");
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
		if(params.parameterExists("filterReadsLongerThanFragment")){
			_longerThanFragmentFilter.filter("Longer than fragment");
			logfile->list("Reads longer than fragment size: filter out. (parameter 'filterReadsLongerThanFragment')");
		} else {
			_longerThanFragmentFilter.keep();
			logfile->list("Reads longer than fragment size: keep. (use 'filterReadsLongerThanFragment' to filter out)");
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
			_firstMateFilter.filter("Second mate");
			_secondMateFilter.keep();
			logfile->list("Mate: keep only first. (parameter 'keepOnlyFirst')");
		}
		else if(params.parameterExists("keepOnlySecond")){
			_firstMateFilter.keep();
			_secondMateFilter.filter("First mate");
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
			logfile->list("Blacklist: keep all. (use 'blacklist' to provide a list and filter specific reads)");
		}

		//Mapping quality filter
		if(params.parameterExists("filterMQ")){
			TNumericRange<uint8_t> Range;
			params.fillParameter("filterMQ", Range);

			_mappingQualityFilter.filter(Range, "Mapping quality outside " + Range.rangeString());
			logfile->list("Mapping quality: restrict to range " + _mappingQualityFilter.rangeString() + ". (parameter 'filterMQ')");
		} else {
			_mappingQualityFilter.keep();
			logfile->list("Mapping quality: keep all. (use 'filterMQ' to limit)");
		}

		//Read length filter
		if(params.parameterExists("filterReadLength")){
			TNumericRange<uint32_t> Range;
			params.fillParameter("filterReadLength", Range);

			_readLengthFilter.filter(Range, "Read length outside " + Range.rangeString());
			logfile->list("Read length: restrict to range " + _readLengthFilter.rangeString() + ". (parameter 'filterReadLength')");
		} else {
			_readLengthFilter.keep();
			logfile->list("Read length: keep all. (use 'filterReadLength' to limit)");
		}


		//Fragment length filter
		if(params.parameterExists("filterFragmentLength")){
			TNumericRange<uint32_t> Range;
			params.fillParameter("filterFragmentLength", Range);

			_fragmentLengthFilter.filter(Range, "Fragment length outside " + Range.rangeString());
			logfile->list("Fragment length: restrict to range " + _fragmentLengthFilter.rangeString() + ". (parameter 'filterFragmentLength')");
		} else {
			_fragmentLengthFilter.keep();
			logfile->list("Fragment length: keep all. (use 'filterFragmentLength' to limit)");
		}
	}
	logfile->endIndent();

	//log filtered reads?
	openBamLog(params, logfile);
};

void TBamFile::curFilterOut(){
	_externalFilter.filterOut(_curBamAlignment.Name, _curBamAlignment.IsReverseStrand(), _curReadGroupID);
};

void TBamFile::filterOut(const std::string & alignmentName, const bool & isReverseStrand, const uint16_t readGroup){
	_externalFilter.filterOut(alignmentName, isReverseStrand, readGroup);
};

void TBamFile::setExternalFilterReason(const std::string reason){
	_externalFilter.setReason(reason);
};

void TBamFile::openBamLog(TParameters & params, TLog* logfile){
	if(params.parameterExists("bamLog") && !_bamLog){
		std::string logFilename = params.getParameter<std::string>("bamLog");
		if(logFilename.empty()){
			logFilename = _filename;
			logFilename = coretools::str::extractBeforeLast(logFilename, ".");
			logFilename += ".bamlog.txt.gz";
		}
		logfile->list("Will write all filtered out reads to '" + logFilename + "'.");
		_bamLog = std::make_shared<TBamFileLog>(logFilename);

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
	if(_bamLog){
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

void TBamFile::_fillChromosomes(genometools::TChromosomes & Chromosomes){
	if(_bamHeader.Sequences.Size() < 1){
		throw "No chromosomes present in BAM header!";
	}

	//make sure object is empty
	Chromosomes.clear();

	//copy from BamHeader
	for(BamTools::SamSequenceIterator chrIt=_bamHeader.Sequences.Begin(); chrIt!=_bamHeader.Sequences.End(); ++chrIt){
		Chromosomes.appendChromosome(chrIt->Name, coretools::str::convertString<uint64_t>(chrIt->Length));
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
	_fileSize = _bamReader.Tell();
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
	if(!_mappedLengthFilter.pass(_curCigar.lengthMapped(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)){
		if(!_allowTooLongReads){
			throw "The mapping length of alignment '" +  _curBamAlignment.Name + "' is beyond the range " + _mappedLengthFilter.rangeString() + "!\n"
			     + "You see this error because " + coretools::__GLOBAL_APPLICATION_NAME__ + " was run with default mapping length filters. Either set your filters using 'filterMappingLength' or add 'allowTooLongReads' to ignore this error.";
		} else {
			_QCFiltersPassed = false;
		}
	} else if(_keepAll){
		//keep all?
		_QCFiltersPassed = true;
	} else {
		//apply regular filters
		_QCFiltersPassed =  _duplicateFilter.pass(!_curBamAlignment.IsDuplicate(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _softClippedFilter.pass(_curCigar.lengthSoftClipped() == 0, _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _improperPairsFilter.pass(!_curBamAlignment.IsPaired() || _curBamAlignment.IsProperPair(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _unmappedFilter.pass(_curBamAlignment.IsMapped(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _failedQCFilter.pass(!_curBamAlignment.IsFailedQC(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _secondaryFilter.pass(_curBamAlignment.IsPrimaryAlignment(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _supplementaryFilter.pass(!_curBamAlignment.IsSupplementary(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _readGroupFilter.pass(_readGroups.readGroupInUse(_curReadGroupID), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _fwdStrandFilter.pass(_curBamAlignment.IsReverseStrand(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _revStrandFilter.pass(!_curBamAlignment.IsReverseStrand(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _firstMateFilter.pass(_curBamAlignment.IsFirstMate(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _secondMateFilter.pass(_curBamAlignment.IsSecondMate(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _mappingQualityFilter.pass(_curBamAlignment.MapQuality, _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _blacklistFilter.pass(!_blacklist.isInBlacklist(_curBamAlignment.Name), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
						 && _readLengthFilter.pass(_curCigar.lengthRead(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID);

		//fragment length
		if(_QCFiltersPassed){
			_QCFiltersPassed = _fragmentLengthFilter.pass(curFragmentLength(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID)
				&& _longerThanFragmentFilter.pass(!_curBamAlignment.IsProperPair() || abs(_curBamAlignment.InsertSize) >= static_cast<int32_t>(_curCigar.lengthAligned()), _curBamAlignment.Name, _curBamAlignment.IsSecondMate(), _curReadGroupID);
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
					throw "Chromosome with refID " + coretools::str::toString(_curBamAlignment.RefID) + " is missing from BAM header!";
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
		throw "BAM file must be sorted by position! Alignment '" + _curBamAlignment.Name + "' is at position " + std::to_string(_curBamAlignment.Position) + ", which is before the position of the previous alignment (" + std::to_string(_previousAlignmentPosition.position()) + ")";
	}

	//store current read group ID
	std::string readGroup;
	_curBamAlignment.GetTag("RG", readGroup);
	_curReadGroupID = _readGroups.getId(readGroup);

	//also update counter per readgroup
	_numAlignmentReadPerReadGroup.add(_curReadGroupID);

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

bool TBamFile::jump(const genometools::TGenomePosition Position){
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

/*
uint16_t TBamFile::curUsableAlignedLength(TQualityFilter & qualFilter) const{
	constexpr char N = genometools::base2char(genometools::Base::N);
	uint16_t counter = 0;
	for(size_t d=0; d<_curBamAlignment.AlignedBases.length(); ++d){
		if(_curBamAlignment.AlignedBases.at(d) != N && qualFilter.pass( BaseQuality(_curBamAlignment.AlignedQualities.at(d)))){
			++counter;
		}
	}
	return counter;
	};*/

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
void TBamFile::printSummaryNoEndIndent(std::string &outputName){
	_logfile->startIndent("Summary of parsed reads from BAM file '" + _filename + "':");
	_logfile->list("Total number of reads read: " + coretools::str::toString(_numAlignmentRead));
	_logfile->list("Reads that passed filters: " + coretools::str::toString(_numAlignmentsPassedQC) + " (" + coretools::str::toPercentString(_numAlignmentsPassedQC, _numAlignmentRead, 3) + "%)");
	uint64_t numFiltered = _numAlignmentRead - _numAlignmentsPassedQC;
	_logfile->list("Reads that were filtered out: " + coretools::str::toString(numFiltered) + " (" + coretools::str::toPercentString(numFiltered, _numAlignmentRead, 3) + "%)");

	if(numFiltered > 0){
		_logfile->addIndent();

				//write counts of filtered reads for each read group to _filterSummary.txt file
				std::string filename = outputName + "_filterSummary.txt";
				coretools::instances::logfile().listFlush("Writing general filter counts to '" + filename + "' ...");

				//creating header
				std::vector<std::string> header;
				header.push_back("readGroup");
				_duplicateFilter.fillHeader(header);
				_softClippedFilter.fillHeader(header);
				_improperPairsFilter.fillHeader(header);
				_unmappedFilter.fillHeader(header);
				_failedQCFilter.fillHeader(header);
				_secondaryFilter.fillHeader(header);
				_supplementaryFilter.fillHeader(header);
				_longerThanFragmentFilter.fillHeader(header);
				_readGroupFilter.fillHeader(header);
				_fwdStrandFilter.fillHeader(header);
				_revStrandFilter.fillHeader(header);
				_firstMateFilter.fillHeader(header);
				_secondMateFilter.fillHeader(header);
				_blacklistFilter.fillHeader(header);
				_mappingQualityFilter.fillHeader(header);
				_fragmentLengthFilter.fillHeader(header);
				_externalFilter.fillHeader(header);
				_readLengthFilter.fillHeader(header);
				_mappedLengthFilter.fillHeader(header);
				coretools::TOutputFile out(filename, header, "\t");

				out << "allReadGroups";
				_duplicateFilter.printCombinedCounts(out);
				_softClippedFilter.printCombinedCounts(out);
				_improperPairsFilter.printCombinedCounts(out);
				_unmappedFilter.printCombinedCounts(out);
				_failedQCFilter.printCombinedCounts(out);
				_secondaryFilter.printCombinedCounts(out);
				_supplementaryFilter.printCombinedCounts(out);
				_longerThanFragmentFilter.printCombinedCounts(out);
				_readGroupFilter.printCombinedCounts(out);
				_fwdStrandFilter.printCombinedCounts(out);
				_revStrandFilter.printCombinedCounts(out);
				_firstMateFilter.printCombinedCounts(out);
				_secondMateFilter.printCombinedCounts(out);
				_blacklistFilter.printCombinedCounts(out);
				_mappingQualityFilter.printCombinedCounts(out);
				_fragmentLengthFilter.printCombinedCounts(out);
				_externalFilter.printCombinedCounts(out);
				_readLengthFilter.printCombinedCounts(out);
				_mappedLengthFilter.printCombinedCounts(out);
				out.endln();


				//writes numbers of removed reads for each applied filter per read group, also lists filters if no reads were removed
				for (uint16_t it = 0; it < _readGroups.size(); it++){
					out << _readGroups.getName(it);
					_duplicateFilter.printCounts(out, it);
					_softClippedFilter.printCounts(out, it);
					_improperPairsFilter.printCounts(out, it);
					_unmappedFilter.printCounts(out, it);
					_failedQCFilter.printCounts(out, it);
					_secondaryFilter.printCounts(out, it);
					_supplementaryFilter.printCounts(out, it);
					_longerThanFragmentFilter.printCounts(out, it);
					_readGroupFilter.printCounts(out, it);
					_fwdStrandFilter.printCounts(out, it);
					_revStrandFilter.printCounts(out, it);
					_firstMateFilter.printCounts(out, it);
					_secondMateFilter.printCounts(out, it);
					_blacklistFilter.printCounts(out, it);
					_mappingQualityFilter.printCounts(out, it);
					_fragmentLengthFilter.printCounts(out, it);
					_externalFilter.printCounts(out, it);
					_readLengthFilter.printCounts(out, it);
					_mappedLengthFilter.printCounts(out, it);
					out.endln();
				}

				out.close();
				coretools::instances::logfile().done();


		//print counts of filtered reads for each read group to terminal, doesn't list filters if no reads were removed
		for (uint16_t it = 0; it < _readGroups.size(); it++){
			_logfile->newLine();
			_logfile->list("Number of reads filtered from read group: '" + coretools::str::toString(_readGroups.getName(it)) + "'");
			_logfile->addIndent();
			_duplicateFilter.summary(_logfile, numFiltered, it);
			_softClippedFilter.summary(_logfile, numFiltered, it);
			_improperPairsFilter.summary(_logfile, numFiltered, it);
			_unmappedFilter.summary(_logfile, numFiltered, it);
			_failedQCFilter.summary(_logfile, numFiltered, it);
			_secondaryFilter.summary(_logfile, numFiltered, it);
			_supplementaryFilter.summary(_logfile, numFiltered, it);
			_longerThanFragmentFilter.summary(_logfile, numFiltered, it);
			_readGroupFilter.summary(_logfile, numFiltered, it);
			_fwdStrandFilter.summary(_logfile, numFiltered, it);
			_revStrandFilter.summary(_logfile, numFiltered, it);
			_firstMateFilter.summary(_logfile, numFiltered, it);
			_secondMateFilter.summary(_logfile, numFiltered, it);
			_blacklistFilter.summary(_logfile, numFiltered, it);
			_mappingQualityFilter.summary(_logfile, numFiltered, it);
			_fragmentLengthFilter.summary(_logfile, numFiltered, it);
			_externalFilter.summary(_logfile, numFiltered, it);
			_readLengthFilter.summary(_logfile, numFiltered, it);
			_mappedLengthFilter.summary(_logfile, numFiltered, it);
			_logfile->endIndent();
		}

		_logfile->endIndent();
	}
};

void TBamFile::printSummary(std::string &outputName){
	printSummaryNoEndIndent(outputName);
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
		_logfile->list("Parsed " + _millionReadsRead() + " million reads (est. " + coretools::str::toStringWithPrecision(positionInFile() * 100, 2) + "%) in " + _timer.formattedTime());
		_lastProgressPrinted = _numAlignmentRead;
	}
};

void TBamFile::printEndWithSummary(std::string &outputName){
	printEndNoEndIndent();
	_logfile->endIndent();
	printSummary(outputName);
};

void TBamFile::printEndNoEndIndent(){
	_logfile->list("Reached end of BAM file in " + _timer.formattedTime() + ':');
	_logfile->conclude("Parsed a total of " + _millionReadsRead() + " million reads in " + _timer.formattedTime() + '.');
};

//------------------------------------------------
// TQualityAdjusterForWriting
//------------------------------------------------
TQualityAdjusterForWriting::TQualityAdjusterForWriting(){
	_adjust = false;
	_binIllumina = false;
	_limitRange = false;
};

TQualityAdjusterForWriting::TQualityAdjusterForWriting(coretools::TParameters & params, coretools::TLog* logfile){
	initialize(params, logfile);
}

void TQualityAdjusterForWriting::initialize(coretools::TParameters & params, coretools::TLog* logfile){
	if(params.parameterExists("outQual")){
		TNumericRange<uint8_t> qualRange;
		params.fillParameter("outQual",  qualRange);
		limitRange(qualRange);

		logfile->list("Will print qualities truncated to ", rangeString(), " (parameter 'outQual')");


		if(qualRange.max() > BaseQuality::max().get()){
			logfile->warning("Truncated quality range to BAM limits!");
		}

	} else {
		logfile->list("Will use the full range of quality scores when writing alignments. (use 'outQual' to constrain).");
	}

	//quality binning
	if(params.parameterExists("writeBinnedQualities")){
		logfile->list("Will write Illumina-binned quality scores. (parameter 'writeBinnedQualities')");

		if(adjusts()){
			logfile->warning("If both 'outQual' and 'writeBinnedQualities' are given, qualities will be truncated first, then binned, and may thus fall outside the requested range.");
		}
		binQualitiesIllumina();
	} else {
		logfile->list("Will write raw quality scores. (use 'writeBinnedQualities' to bin)");
	}
};

void TQualityAdjusterForWriting::binQualitiesIllumina(){
	_binIllumina = true;
	_adjust = true;
};

void TQualityAdjusterForWriting::limitRange(const BaseQuality & min, const BaseQuality & max){
	_minQual = min;
	_maxQual = max;
	_adjust = true;
};

void TQualityAdjusterForWriting::limitRange(const TNumericRange<uint8_t> & Range){
	if(Range.minIncluded()){
		_minQual = Range.min();
	} else {
		_minQual = Range.min() + 1;
	}
	if(Range.maxIncluded()){
		_maxQual = Range.min();
	} else {
		_maxQual = Range.min() - 1;
	}
	_adjust = true;
};

std::string TQualityAdjusterForWriting::rangeString(){
	return "[" + coretools::str::toString(genometools::PhredIntProbability(_minQual)) + "," + coretools::str::toString(genometools::PhredIntProbability(_maxQual)) + "]";
};

char TQualityAdjusterForWriting::_adjustOneQuality(BaseQuality qual) const {
	if(qual < _minQual){
		qual = _minQual;
	} else if (qual > _maxQual){
		qual = _maxQual;
	}
	if(_binIllumina){
		qual.makeIllumina();
	}

	return qual.get();
};

void TQualityAdjusterForWriting::adjustQualities(std::string & qualities) const {
	if(_adjust){
		for(auto& q : qualities){
			q = _adjustOneQuality(BaseQuality(q));
		}
	}
};

//-----------------------------------------------------
//TOutputBamFile
//----------------------------------------------------
TOutputBamFile::TOutputBamFile(){
	_openForWriting = false;
	_readGroups = nullptr;
};

TOutputBamFile::TOutputBamFile(const std::string Filename, const TBamFile & Original){
	_openForWriting = false;
	open(Filename, Original.samHeader(), Original.chromosomes(), Original.readGroups());
};

TOutputBamFile::~TOutputBamFile(){
	closeNoIndex();
};

void TOutputBamFile::open(const std::string Filename, const TSamHeader & Header, const genometools::TChromosomes & Chromosomes, const TReadGroups & ReadGroups){
	closeNoIndex();

	_outputFilename = Filename;
	_readGroups = &ReadGroups;

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

void TOutputBamFile::open(const std::string Filename, const TBamFile & Original){
	open(Filename, Original.samHeader(), Original.chromosomes(), Original.readGroups());
};

void TOutputBamFile::open(TParameters & params, TLog* logfile, const std::string Filename, const TSamHeader & Header, const genometools::TChromosomes & Chromosomes, const TReadGroups & ReadGroups){
	logfile->startIndent("Writing alignments to new BAM to file '" + Filename + "':");

	//first open bam file
	open(Filename, Header, Chromosomes, ReadGroups);

	//read output settings
	setQualityAdjusterForWriting(params, logfile);
};

void TOutputBamFile::setQualityAdjusterForWriting(coretools::TParameters & params, coretools::TLog* logfile){
	_qualityAdjuster.initialize(params, logfile);
};

void TOutputBamFile::setQualityAdjusterForWriting(const TQualityAdjusterForWriting & QualityAdjuster){
	_qualityAdjuster = QualityAdjuster;
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
	_qualityAdjuster.adjustQualities(alignment.Qualities);

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
	_tmpBamAlignment.QueryBases = alignment.sequence();
	_tmpBamAlignment.Qualities = alignment.qualities();

	//add read group information
	_tmpBamAlignment.AddTag("RG", "Z", _readGroups->getName(alignment.readGroupId()));

	//and now write
	if(!_bamWriter.SaveAlignment(_tmpBamAlignment)){
		throw "Read '" + _tmpBamAlignment.Name + "' could not be written!";
	}
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

