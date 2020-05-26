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
	_maxReadLength = 200;
	_fileSize = 0;
	_numAlignmentRead = 0;
	_numAlignmentsPassedQC = 0;
	_limitNumReads = false;
	_maxNumReadsToRead = 0;

	//current alignment
	_curReadGroupID = 0;
	_previousAlignmentPos = 0;
	_previousAlignmentChr = -1;
	_chrChanged = false;

	//writing
	_openForWriting = false;

	//set filters to default
	_QCFiltersPassed = false;
	_maxReadLength = 65536;
	_keepAll = false;

	//blacklist
	_updateLog = false;
	_bamLog = nullptr;

	//progress reporting
	_logfile = nullptr;
	_progressFrequency = 100000;
	_lastProgressPrinted = 0;
};

void TBamFile::setFiltersAndLimits(TParameters & params, TLog* logfile){
	_limitChromosomes(params, logfile);
	_limitReadGroups(params, logfile);
	_setFilters(params, logfile);
};

void TBamFile::_limitReadGroups(TParameters & params, TLog* logfile){
	if(params.parameterExists("readGroup")){
		readGroups.filterReadGroups(params.getParameterString("readGroup"));
		logfile->startIndent("Will limit analysis to the following read groups:");
		readGroups.printReadgroupsInUse(logfile);
		logfile->endIndent();
		_readGroupFilter.filter("Read group not in use");
	} else {
		_readGroupFilter.keep();
	}
};

void TBamFile::limitReadLength(const int MaxReadLength){
	//set max read length, must be >=1 but smaller than 65536 (uint16_t)
	if(MaxReadLength < 1)
		throw "Max read length must be at least 1 bp!";
	if(MaxReadLength > 65536)
		throw "Atlas currently only supports read length up to 65536 bp. Contact us if you have longer reads / fragments";
	_maxReadLength = MaxReadLength;
};

void TBamFile::_setFilters(TParameters & params, TLog* logfile){
	//max read length
	int MaxReadLength = params.getParameterIntWithDefault("maxReadLength", 200);
	logfile->list("Expect no read to be longer than " + toString(MaxReadLength) + ". (parameter 'maxReadLength')");
	limitReadLength(MaxReadLength);

	//number of reads
	if(params.parameterExists("limitReads")){
		_maxNumReadsToRead = params.getParameterInt("limitReads");
		logfile->list("Will limit the analysis to the first " + toString(_maxNumReadsToRead) + " reads in the BAM file.");
		_limitNumReads = true;
	}

	//alignment filters
	logfile->startIndent("Will use the following filters on reads:");
	if(params.parameterExists("keepAllReads")){
		_keepAll = true;
		logfile->list("Will keep all reads. (parameter 'keepDuplicates', overrules any other QC filter)");
	} else {
		_keepAll = false;
		//duplicates
		if(params.parameterExists("keepDuplicates")){
			_duplicateFilter.keep();
			logfile->list("Will keep duplicate reads. (parameter 'keepDuplicates')");
		} else {
			_duplicateFilter.filter("Duplicate");
			logfile->list("Will filter out duplicate reads. (use 'keepDuplicates' to keep)");
		}

		//soft clips
		if(params.parameterExists("filterSoftClips")){
			_softClippedFilter.filter("Soft clipped");
			logfile->list("Will filter out soft clipped reads. (parameter 'filterSoftClips')");
		} else {
			_softClippedFilter.keep();
			logfile->list("Will keep soft clipped reads. (use 'filterSoftClips' to filter out)");
		}

		//improper pairs
		if(params.parameterExists("keepImproperPairs")){
			_improperPairsFilter.keep();
			logfile->list("Will keep improper pairs. (parameter 'keepImproperPairs')");
		} else {
			_improperPairsFilter.filter("Improper pair");
			logfile->list("Will filter out improper pairs. (use 'keepImproperPairs' to keep)");
		}

		//unmapped reads
		if(params.parameterExists("keepUnmappedReads")){
			_unmappedFilter.keep();
			logfile->list("Will keep unmapped reads. (parameter 'keepUnmappedReads')");
		} else {
			_unmappedFilter.filter("Unmapped");
			logfile->list("Will filter out unmapped reads. (use 'keepUnmappedReads' to keep)");
		}

		//failed QC
		if(params.parameterExists("keepFailedQC")){
			_failedQCFilter.keep();
			logfile->list("Will keep reads that failed QC. (parameter 'keepFailedQC')");
		} else {
			_failedQCFilter.filter("Failed QC");
			logfile->list("Will filter out reads that failed QC. (use 'keepFailedQC' to keep)");
		}

		//secondary reads
		if(params.parameterExists("keepSecondaryReads")){
			_secondaryFilter.keep();
			logfile->list("Will keep secondary reads. (parameter 'keepSecondaryReads')");
		} else {
			_secondaryFilter.filter("Secondary alignment");
			logfile->list("Will filter out secondary reads. (use 'keepSecondaryReads' to keep)");
		}

		//supplementary reads
		if(params.parameterExists("keepSupplementaryReads")){
			_supplementaryFilter.keep();
			logfile->list("Will keep supplementary reads. (parameter 'keepSupplementaryReads')");
		} else {
			_supplementaryFilter.filter("Supplementary alignment");
			logfile->list("Will filter out supplementary reads. (use 'keepSupplementaryReads' to keep)");
		}

		//fragment length
		if(params.parameterExists("keepReadsLongerThanFragment")){
			_longerThanFragmentFilter.keep();
			logfile->list("Will keep reads that are longer than the fragment size. (parameter 'keepReadsLongerThanFragment')");
		} else {
			_longerThanFragmentFilter.filter("Longer than fragment");
			logfile->list("Will filter out reads that are longer than the fragment size. (use 'keepReadsLongerThanFragment' to keep)");
		}

		//strand
		if(params.parameterExists("keepOnlyFwd")){
			_fwdStrandFilter.keep();
			_revStrandFilter.filter("Reverse strand");
			logfile->list("Will keep only forward mapping reads. (parameter 'keepOnlyFwd')");
		}
		else if(params.parameterExists("keepOnlyRev")){
			_fwdStrandFilter.filter("Forward strand");
			_revStrandFilter.keep();
			logfile->list("Will keep only reverse mapping reads. (parameter 'keepOnlyRev')");
		} else {
			_fwdStrandFilter.keep();
			_revStrandFilter.keep();
			logfile->list("Will keep forward and reverse mapping reads. (use 'keepOnlyFwd' or 'keepOnlyRev' to limit)");
		}

		//mate
		if(params.parameterExists("keepOnlyFirst")){
			_firstMateFilter.keep();
			_secondMateFilter.filter("Second mate");
			logfile->list("Will keep only the first mates. (parameter 'keepOnlyFirst')");
		}
		else if(params.parameterExists("keepOnlySecond")){
			_firstMateFilter.filter("First mate");
			_secondMateFilter.keep();
			logfile->list("Will keep only the second mates. (parameter 'keepOnlySecond')");
		} else {
			_firstMateFilter.keep();
			_secondMateFilter.keep();
			logfile->list("Will keep first and second mates. (use 'keepOnlyFirst' or 'keepOnlySecond' to limit)");
		}

		//blacklist
		if(params.parameterExists("blacklist")){
			std::string blaclistFilename = params.getParameterString("blacklist");
			logfile->list("Will filter out reads present in the file '" + blaclistFilename + "'. (parameter 'blacklist')");
			_blacklist.addFromFile(blaclistFilename);
			_blacklistFilter.filter("was in provided blacklist");
		} else {
			_blacklistFilter.keep();
		}

		//Mapping quality filter
		if(params.parameterExists("minMQ") || params.parameterExists("maxMQ")){
			int MinMQ = params.getParameterInt("minMQ", 0);
			int MaxMQ = params.getParameterInt("maxMQ", 254);

			if(MinMQ < 0 || MinMQ > 254)
				throw "minMQ '" + toString(MinMQ) + "' is outside the accepted range [0,254]!";
			if(MaxMQ < 0 || MaxMQ > 254)
				throw "maxMQ '" + toString(MaxMQ) + "' is outside the accepted range [0,254]!";
			if(MaxMQ < MinMQ)
				throw "maxMQ must be <= minMQ";

			_mappingQualityFilter.filter(MinMQ, MaxMQ, "Mapping quality outside [" + toString(MinMQ) + ", " + toString(MaxMQ) + "]");
			logfile->list("Will filter out reads with mapping quality outside the range [" + toString(MinMQ) + ", " + toString(MaxMQ) + "]. (parameters 'minMQ', 'maxMQ')");
		} else {
			_mappingQualityFilter.keep();
			logfile->list("Will keep reads regardless of their mapping quality. (use 'minMQ' and 'maxMQ' to limit)");
		}

		//Fragment length filter
		if(params.parameterExists("minFragmentLength") || params.parameterExists("maxFragmentLength")){
			int MinFragmentLength = params.getParameterIntWithDefault("minFragmentLength", 0);
			int MaxFragmentLength = params.getParameterIntWithDefault("maxFragmentLength", 1000);
			if(MinFragmentLength < 0)
				throw "minFragmentLength '" + toString(MinFragmentLength) + "' must be >0 0!";
			if(MinFragmentLength < 0)
				throw "minFragmentLength '" + toString(MinFragmentLength) + "' must be >0 0!";
			if(MinFragmentLength > MaxFragmentLength)
				throw "minFragmentLength must be <= maxFragmentLength!";

			_fragmentLengthfilter.filter(MinFragmentLength, MaxFragmentLength, "Fragment length outside [" + toString(MinFragmentLength) + ", " + toString(MaxFragmentLength) + "]");
			logfile->list("Will filter out reads with fragment length outside the range [" + toString(MinFragmentLength) + ", " + toString(MaxFragmentLength) + "]. (parameters 'minFragmentLength', 'maxFragmentLength')");
		} else {
			_fragmentLengthfilter.keep();
			logfile->list("Will keep reads reagrless of their fragment length. (use 'minFragmentLength', 'maxFragmentLength' to limit)");
		}
	}
	logfile->endIndent();

	//log filtered reads?
	openBamLog(params, logfile);
};

void TBamFile::openBamLog(TParameters & params, TLog* logfile){
	if(params.parameterExists("bamLog") && !_updateLog){
		std::string logFilename = params.getParameterString("bamLog");
		if(logFilename.empty()){
			logFilename = _filename;
			logFilename = extractBeforeLast(logFilename, ".");
			logFilename += ".bamlog.txt.gz";
		}
		logfile->list("Will all filtered out reads to '" + logFilename + "'.");
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
		_fragmentLengthfilter.setLog(_bamLog);
	}
};

void TBamFile::_fillChromosomes(TChromosomes & chromosomes){
	if(_bamHeader.Sequences.Size() < 1){
		throw "No chromosomes present in BAM header!";
	}

	//make sure object is empty
	chromosomes.clear();

	//copy from BamHeader
	for(BamTools::SamSequenceIterator chrIt=_bamHeader.Sequences.Begin(); chrIt!=_bamHeader.Sequences.End(); ++chrIt){
		chromosomes.appendChromosome(chrIt->Name, stringToLong(chrIt->Length));
	}
};

//--------------------------------------------------------
// Functions for reading
//--------------------------------------------------------
void TBamFile::open(const std::string filename, const bool indexNotRequired){
	//open BAM file
	if (!_bamReader.Open(filename))
		throw "Failed to open BAM file '" + filename + "'!";

	//load index file
	if(!_bamReader.LocateIndex() && !indexNotRequired)
		throw "No index file found for BAM file '" + filename + "'!";

	//initialize bam stuff
	_bamHeader = _bamReader.GetHeader();

	//initialize read groups
	readGroups.fill(_bamHeader);

	//initialize chromosomes
	_fillChromosomes(chromosomes);

	//parse CIGAR
	_curCigar.clear();
	for(auto& it : _curBamAlignment.CigarData){
		_curCigar.add(it.Type, it.Length);
	}

	//get file size
	_bamReader.Jump(chromosomes.size() - 1, 0);
	BamTools::BamAlignment bamAlignment;
	_bamReader.GetNextAlignment(bamAlignment);
	_fileSize = _bamReader.tell();
	_bamReader.Rewind();
};

void TBamFile::_applyFilters(){
	//keep all?
	if(_keepAll){
		_QCFiltersPassed = true;
	} else {
		_QCFiltersPassed =  _duplicateFilter.pass(_curBamAlignment.IsDuplicate(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _softClippedFilter.pass(_curCigar.lengthSoftClipped() > 0, _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _improperPairsFilter.pass(_curBamAlignment.IsPaired() && !_curBamAlignment.IsProperPair(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _unmappedFilter.pass(!_curBamAlignment.IsMapped(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _failedQCFilter.pass(_curBamAlignment.IsFailedQC(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _secondaryFilter.pass(!_curBamAlignment.IsPrimaryAlignment(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _supplementaryFilter.pass(_curBamAlignment.IsSupplementary(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _readGroupFilter.pass(!readGroups.readGroupInUse(_curReadGroupID), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _fwdStrandFilter.pass(!_curBamAlignment.IsReverseStrand(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _revStrandFilter.pass(_curBamAlignment.IsReverseStrand(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _firstMateFilter.pass(_curBamAlignment.IsFirstMate(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _secondMateFilter.pass(_curBamAlignment.IsSecondMate(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _mappingQualityFilter.pass(_curBamAlignment.MapQuality, _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
						 && _blacklistFilter.pass(_blacklist.isInBlacklist(_curBamAlignment.Name), curBamAlignment.Name, _curBamAlignment.IsSecondMate());

		//fragment length
		if(_QCFiltersPassed){
			uint32_t fragmentLength;
			if(_curBamAlignment.IsProperPair()){
				fragmentLength = abs(_curBamAlignment.InsertSize) + _curCigar.lengthInserted() - _curCigar.lengthDeleted();
			} else {
				fragmentLength = _curCigar.lengthSequenced();
			}
			_QCFiltersPassed = _fragmentLengthfilter.pass(fragmentLength, _curBamAlignment.Name, _curBamAlignment.IsSecondMate())
					        && _longerThanFragmentFilter.pass(_curBamAlignment.InsertSize < _curCigar.lengthAligned(), _curBamAlignment.Name, _curBamAlignment.IsSecondMate());
		}
	}

	//update counter
	if(_QCFiltersPassed){
		++_numAlignmentsPassedQC;
	}
};

bool TBamFile::readNextAlignment(){
	if(_limitNumReads && _numAlignmentRead >=_maxNumReadsToRead){
		return false;
	}

	if(!_bamReader.GetNextAlignment(_curBamAlignment)){
		return false;
	}

	//check if chromosome changed
	if(_curBamAlignment.RefID != chromosomes.curRefID()){
		_previousAlignmentPos = 0;
		_previousAlignmentChr = _curBamAlignment.RefID;

		//advance chromosomes
		chromosomes.jumpTo(_curBamAlignment.RefID);

		//if not in use: jump to next in use
		if(!chromosomes.curInUse()){
			while(!chromosomes.curInUse()){
				if(!chromosomes.next()){
					return false;
				}
			}

			//jump reader and read first alignment
			jump(chromosomes.curRefID(), 0);
			if(!_bamReader.GetNextAlignment(_curBamAlignment)){
				return false;
			}
		}
		_chrChanged = true;
	} else {
		_chrChanged = false;
	}

	++_numAlignmentRead;

	//check if BAM file is sorted
	if(_curBamAlignment.Position < _previousAlignmentPos)
		throw "BAM file must be sorted by position! Alignment '" + _curBamAlignment.Name + "' is at position " + toString(_curBamAlignment.Position) + ", which is before the position of the previous alignment (" + toString(_previousAlignmentPos) + ")";
	_previousAlignmentPos = _curBamAlignment.Position;

	//check read length
	if(_curBamAlignment.AlignedBases.size() > _maxReadLength)
		throw "Alignment '" +  _curBamAlignment.Name + "' is longer than the max read length " + toString(_maxReadLength) + "! Please change max read length to parse this data.";

	//store current read group ID
	std::string readGroup;
	_curBamAlignment.GetTag("RG", readGroup);
	_curReadGroupID = readGroups.find(readGroup);

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

void TBamFile::fill(TAlignment & alignment){
	alignment.fill(_curBamAlignment.Name,
				   BAM::TSamFlags(_curBamAlignment.AlignmentFlag),
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

bool TBamFile::jump(int chr, int position){
	_previousAlignmentPos = -1;
	return _bamReader.Jump(chr, position);
};

void TBamFile::rewind(){
	_bamReader.Rewind();
	_numAlignmentRead = 0;
	_numAlignmentsPassedQC = 0;
	_previousAlignmentPos = -1;
};

//--------------------------------------------------------
// Functions for writing
//--------------------------------------------------------
void TBamFile::openOutput(std::string filename){
	BamTools::RefVector references = _bamReader.GetReferenceData();
	if(!bamWriter.Open(filename, _bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";
	_openForWriting = true;
};

void TBamFile::writeCurAlignment(){
	if(!_openForWriting){
		throw "BAM writer is not open!";
	}

	// write alignment
	if(!bamWriter.SaveAlignment(_curBamAlignment))
		throw "Read '" + _curBamAlignment.Name + "' could not be written!";
};

void TBamFile::writeAlignment(TAlignment & alignment, const TGenotypeMap & genoMap, const TQualityMap & qualityMap){
	if(!_openForWriting){
		throw "BAM writer is not open!";
	}

	//create bamAlignment and then write
	BamTools::BamAlignment _tmpBamAlignment;

	_tmpBamAlignment.Name = alignment.name;
	_tmpBamAlignment.AlignmentFlag = alignment.flags.asInt();
	_tmpBamAlignment.RefID = alignment.refID;
	_tmpBamAlignment.Position = alignment.position;

	_tmpBamAlignment.MapQuality = alignment.mappingQuality;
	_tmpBamAlignment.MateRefID = alignment.mateRefID;
	_tmpBamAlignment.MatePosition = alignment.matePosition;
	_tmpBamAlignment.InsertSize = alignment.insertSize_TLEN;

	//CIGAR
	for(auto& it : alignment.cigar){
		_tmpBamAlignment.CigarData.emplace_back(it.type, it.length);
	}

	//add sequences and qualities
	_tmpBamAlignment.QueryBases = alignment.getSequence(genoMap, qualityMap);
	_tmpBamAlignment.Qualities = alignment.getQualities(genoMap, qualityMap);

	//add read group information
	_tmpBamAlignment.AddTag("RG", "Z", readGroups.getName(alignment.readGroupID));

	//and now write
	if(!bamWriter.SaveAlignment(_tmpBamAlignment))
		throw "Read '" + _tmpBamAlignment.Name + "' could not be written!";
};

//--------------------------------------------------------
// Getters
//--------------------------------------------------------
int TBamFile::curUsableLength(const int minQual, const int maxQual){
	int counter = _curBamAlignment.Length;
	for(int d=0; d<_curBamAlignment.Length; ++d){
		if(_curBamAlignment.QueryBases.at(d) == N || _curBamAlignment.Qualities.at(d) < minQual || _curBamAlignment.Qualities.at(d) > maxQual){
			counter -= 1;
		}
	}

	return counter;
};

//-----------------------------------------------------
// Reporting
//-----------------------------------------------------
void TBamFile::printSummary(TLog* Logfile){
	Logfile->startIndent("Summary of parsed reads from BAM file '" + _filename + "':")
	Logfile->list("Total number of reads read: " + toString(_numAlignmentRead));
	Logfile->list("Reads that passed filters: " + toString(_numAlignmentsPassedQC));

	Logfile->startIndent(toString(_numAlignmentRead - _numAlignmentsPassedQC) + " reads were filtered out:");

	_duplicateFilter.summary(Logfile);
	_softClippedFilter.summary(Logfile);
	_improperPairsFilter.summary(Logfile);
	_unmappedFilter.summary(Logfile);
	_failedQCFilter.summary(Logfile);
	_secondaryFilter.summary(Logfile);
	_supplementaryFilter.summary(Logfile);
	_longerThanFragmentFilter.summary(Logfile);
	_readGroupFilter.summary(Logfile);
	_fwdStrandFilter.summary(Logfile);
	_revStrandFilter.summary(Logfile);
	_firstMateFilter.summary(Logfile);
	_secondMateFilter.summary(Logfile);
	_blacklistFilter.summary(Logfile);
	_mappingQualityFilter.summary(Logfile);
	_fragmentLengthfilter.summary(Logfile);

	Logfile->endIndent();
	Logfile->endIndent();
};

void TBamFile::startProgressReporting(uint32_t Frequency, TLog* Logfile){
	_logfile = Logfile;
	_progressFrequency = Frequency;
	_lastProgressPrinted = 0;
	_timer.start();

	_logfile->startIndent("Parsing through BAM file:");
};

void TBamFile::printProgress(){
	if(_numAlignmentRead - _lastProgressPrinted >= _progressFrequency){
		_logfile->list("Parsed " + _millionReadsRead() + " million reads (est. " + to_string_with_precision(positionInFile() * 100, 2) + "%) in " + _timer.minutes() + " min.");
		_lastProgressPrinted = _numAlignmentRead;
	}
};

void TBamFile::printEnd(){
	printEndNoEndIndent();
	_logfile->endIndent();
	printSummary(_logfile);
};

void TBamFile::printEndNoEndIndent(){
	_logfile->list("Reached end of BAM file in " + _timer.minutes() + " min.");
	_logfile->conclude("Parsed a total of " + _millionReadsRead() + " million reads in " + _timer.minutes() + " min.");
};

}; //end namespace

