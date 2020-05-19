/*
 * TBamFile.cpp
 *
 *  Created on: May 18, 2020
 *      Author: wegmannd
 */

#include "TBamFile.h"

//-----------------------------------------------------
//TBamFile
//-----------------------------------------------------
TBamFile::TBamFile(){
	_maxReadLength = 200;
	_fileSize = 0;
	_numAlignmentRead = 0;
	_numAlignmentsPassedQC = 0;
	_previousAlignmentPos = -1;
	_previousAlignmentChr = -1;
	_chrChanged = false;
};

void TBamFile::setAlignmentFiltersToKeepAll(){
	_keepAll = true;
	_keepDuplicates = true;
	_filterSoftClips = false;
	_keepImproperPairs = true;
	_keepUnmappedReads = true;
	_keepFailedQC = true;
	_keepSecondary = true;
	_keepSupplementary = true;
	_keepReadsLongerThanInsertSize = true;
};

void TBamFile::setAlignmentFilters(TParameters & params, TLog* logfile){
	//Mapping quality filter
	if(params.parameterExists("minMQ") || params.parameterExists("maxMQ")){
		_setMappingQualityFilters(params.getParameterInt("minMQ"), params.getParameterInt("maxMQ"));
	}

	//Fragment length filter:: CAN I KEEP HERE??? NEED PARSING!!
	if(params.parameterExists("minFragmentLength") || params.parameterExists("maxFragmentLength")){

		int MinFragmentLength = params.getParameterIntWithDefault("minFragmentLength", 0);
		int MaxFragmentLength = params.getParameterIntWithDefault("maxFragmentLength", 255);
		if(MinFragmentLength < 0 || MinFragmentLength > 255)
			throw "minFragmentLength '" + toString(MinFragmentLength) + "' is outside the accepted range [0,255]!";
		if(MinFragmentLength < 0 || MinFragmentLength > 255)
			throw "maxFragmentLength '" + toString(MinFragmentLength) + "' is outside the accepted range [0,255]!";
		if(MinFragmentLength > MaxFragmentLength)
			throw "minFragmentLength must be <= maxFragmentLength!";

		minFragmentLength = MinFragmentLength; maxFragmentLength = MaxFragmentLength;
		logfile->list("Will filter out reads with fragment length outside the range [" + toString(MinFragmentLength) + ", " + toString(MaxFragmentLength) + "]. (parameters 'minFragmentLength', 'maxFragmentLength')");
	}

	//duplicates
	if(params.parameterExists("keepDuplicates")){
		_keepDuplicates = true;
		logfile->list("Will keep duplicate reads. (parameter 'keepDuplicates')");
	}

	//soft clips
	if(params.parameterExists("filterSoftClips")){
		_filterSoftClips = true;
		logfile->list("Will filter out soft clipped reads. (parameter 'filterSoftClips')");
	}

	//improper pairs
	if(params.parameterExists("keepImproperPairs")){
		_keepImproperPairs = true;
		logfile->list("Will keep improper pairs. (parameter 'keepImproperPairs')");
	}

	//unmapped reads
	if(params.parameterExists("keepUnmappedReads")){
		_keepUnmappedReads = true;
		logfile->list("Will keep unmapped reads. (parameter 'keepUnmappedReads')");
	}

	//failed QC
	if(params.parameterExists("keepFailedQC")){
		_keepFailedQC = true;
		logfile->list("Will keep reads that failed QC. (parameter 'keepFailedQC')");
	}

	//secondary reads
	if(params.parameterExists("keepSecondaryReads")){
		_keepSecondary = true;
		logfile->list("Will keep secondary reads. (parameter 'keepSecondaryReads')");
	}

	//supplementary reads
	if(params.parameterExists("keepSupplementaryReads")){
		_keepSupplementary = true;
		logfile->list("Will keep supplementary reads. (parameter 'keepSupplementaryReads')");
	}

	//fragment length
	if(params.parameterExists("keepReadsLongerThanFragment")){
		_keepReadsLongerThanInsertSize = true;
		logfile->list("Will keep reads that are longer than the fragment size. (parameter 'keepReadsLongerThanFragment')");
	}

	//strand
	if(params.parameterExists("keepOnlyFwd")){
		useStrand[1] = false;
		logfile->list("Will keep only forward mapping reads. (parameter 'keepOnlyFwd')");
	}
	else if(params.parameterExists("keepOnlyRev")){
		useStrand[0] = false;
		logfile->list("Will keep only reverse mapping reads. (parameter 'keepOnlyRev')");
	}

	//mate
	if(params.parameterExists("keepOnlyFirst")){
		useMate[1] = false;
		logfile->list("Will keep only the first mates. (parameter 'keepOnlyFirst')");
	}
	else if(params.parameterExists("keepOnlySecond")){
		useMate[0] = false;
		logfile->list("Will keep only the second mates. (parameter 'keepOnlySecond')");
	}
};

void TBamFile::limitReadGroups(std::string readGroupList, TLog* logfile){
	readGroups.filterReadGroups(readGroupList);
	logfile->startIndent("Will limit analysis to the following read groups:");
	readGroups.printReadgroupsInUse(logfile);
	logfile->endIndent();
};

void TBamFile::open(std::string filename, uint32_t maxReadLength, bool indexNotRequired){
	//open BAM file
	if (!_bamReader.Open(filename))
		throw "Failed to open BAM file '" + filename + "'!";
	//load index file
	if(!_bamReader.LocateIndex() && !indexNotRequired)
		throw "No index file found for BAM file '" + filename + "'!";

	//initialize bam stuff
	_bamHeader = _bamReader.GetHeader();

	//set max read length, must be >=1 but smaller than 65536 (uint16_t)
	if(maxReadLength < 1)
		throw "Max read length must be at least 1 bp!";
	if(maxReadLength > 65536)
		throw "Atlas currently only supports read length up to 65536 bp. Contact us if you have longer reads / fragments";
	_maxReadLength = maxReadLength;

	//initialize chromosomes
	chromosomes.readChromosomes(&_bamHeader);

	//initialize read groups
	readGroups.fill(_bamHeader);

	//get file size
	_bamReader.Jump(chromosomes.size() - 1, 0);
	BamTools::BamAlignment bamAlignment;
	_bamReader.GetNextAlignment(bamAlignment);
	_fileSize = _bamReader.tell();
	_bamReader.Rewind();
};

void TBamFile::_applyFilters(){
	std::vector<int> clipSizes;
	std::vector<int> readPositions;
	std::vector<int> genomePositions;

	//standard QC
	_QCFiltersPassed = readGroups.readGroupInUse(_curReadGroupID)
					&& (_keepImproperPairs || (!bamAlignment.IsPaired() || bamAlignment.IsProperPair()))
					&& (_keepUnmappedReads || bamAlignment.IsMapped())
					&& (_keepFailedQC || !bamAlignment.IsFailedQC())
					&& (_keepSecondary || bamAlignment.IsPrimaryAlignment())
					&& (_keepSupplementary || !bamAlignment.IsSupplementary())
					&& chromosomes.inUse(bamAlignment.RefID)
					&& (_keepDuplicates || !bamAlignment.IsDuplicate())
					&& (!_filterSoftClips || !bamAlignment.GetSoftClips(clipSizes, readPositions, genomePositions))
					&& useStrand[bamAlignment.IsReverseStrand()]
					&& useMate[bamAlignment.IsSecondMate()]
					&& (!applyMQFilter || (bamAlignment.MapQuality >= minMQ && bamAlignment.MapQuality <= maxMQ));

	//fragment length
	if(bamAlignment.IsPaired() && abs(bamAlignment.InsertSize) < bamAlignment.AlignedBases.length()){
		_QCFiltersPassed = false;
	}

	if(_QCFiltersPassed){
		++_numAlignmentsPassedQC;
	}
};

bool TBamFile::readNextAlignment(){
	if(!_bamReader.GetNextAlignment(bamAlignment)){
		return false;
	}

	++_numAlignmentRead;

	//check if BAM file is sorted
	if(bamAlignment.Position < _previousAlignmentPos)
		throw "BAM file must be sorted by position! Alignment '" + bamAlignment.Name + "' is at position " + toString(bamAlignment.Position) + ", which is before the position of the previous alignment (" + toString(_previousAlignmentPos) + ")";
	_previousAlignmentPos = bamAlignment.Position;

	//check if chromosome changed
	if(bamAlignment.RefID != _previousAlignmentChr){
		_previousAlignmentPos = 0;
		_previousAlignmentChr = bamAlignment.RefID;
		_chrChanged = true;
	} else {
		_chrChanged = false;
	}

	//check read length
	if(bamAlignment.AlignedBases.size() > _maxReadLength)
		throw "Alignment '" +  bamAlignment.Name + "' is longer than the max read length " + toString(_maxReadLength) + "! Please change max read length to parse this data.";

	//store current read group ID
	std::string readGroup;
	bamAlignment.GetTag("RG", readGroup);
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
	//make sure container is empty
	alignment.clear();

	//add basic info
	alignment.bamAlignment = bamAlignment; //TODO: find way to not rely on bamtools in TAlignment

	//TODO: use functions as for read group?
	alignment.chrNumber = bamAlignment.RefID;
	alignment.position = bamAlignment.Position;
	alignment.alignmentName = bamAlignment.Name;
	alignment.isReverseStrand = bamAlignment.IsReverseStrand();
	alignment.isPaired = bamAlignment.IsPaired();
	alignment.isProperPair = bamAlignment.IsProperPair();
	alignment.mappingQuality = bamAlignment.MapQuality;
	alignment.insertSize = bamAlignment.InsertSize;
	alignment.isSecondMate = bamAlignment.IsSecondMate();
	alignment.matePosition = bamAlignment.MatePosition;

	alignment.fillReadGroupInfo(_curReadGroupID);

	alignment.empty = false;
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

bool TBamFile::curIsLongerThanInsertSize(){
	return bamAlignment.IsPaired() && abs(bamAlignment.InsertSize) < bamAlignment.AlignedBases.length();
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


