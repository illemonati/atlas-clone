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
	_previousAlignmentPos = -1;
	_previousAlignmentChr = -1;
	_chrChanged = false;
};

void TBamFile::open(std::string filename, uint32_t maxReadLength, bool indexNotRequired){
	//open BAM file
	if (!bamReader.Open(filename))
		throw "Failed to open BAM file '" + filename + "'!";
	//load index file
	if(!bamReader.LocateIndex() && !indexNotRequired)
		throw "No index file found for BAM file '" + filename + "'!";

	//initialize bam stuff
	bamHeader = bamReader.GetHeader();

	//set max read length, must be >=1 but smaller than 65536 (uint16_t)
	if(maxReadLength < 1)
		throw "Max read length must be at least 1 bp!";
	if(maxReadLength > 65536)
		throw "Atlas currently only supports read length up to 65536 bp. Contact us if you have longer reads / fragments";
	_maxReadLength = maxReadLength;

	//initialize chromosomes
	chromosomes.readChromosomes(&bamHeader);

	//initialize read groups
	readGroups.fill(bamHeader);

	//get file size
	bamReader.Jump(chromosomes.size() - 1, 0);
	BamTools::BamAlignment bamAlignment;
	bamReader.GetNextAlignment(bamAlignment);
	_fileSize = bamReader.tell();
	bamReader.Rewind();
};

void TBamFile::limitReadGroups(std::string readGroupList, TLog* logfile){
	readGroups.filterReadGroups(readGroupList);
	logfile->startIndent("Will limit analysis to the following read groups:");
	readGroups.printReadgroupsInUse(logfile);
	logfile->endIndent();
};

bool TBamFile::readNextAlignment(){
	if(!bamReader.GetNextAlignment(bamAlignment)){
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

void TBamFile::_applyFilters(){
	std::vector<int> clipSizes;
	std::vector<int> readPositions;
	std::vector<int> genomePositions;

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
};

bool TBamFile::curIsLongerThanInsertSize(){
	return bamAlignment.IsPaired() && abs(bamAlignment.InsertSize) < bamAlignment.AlignedBases.length();
};

bool TBamFile::jump(int chr, int position){
	_previousAlignmentPos = -1;
	return bamReader.Jump(chr, position);
};

void TBamFile::rewind(){
	bamReader.Rewind();
	_numAlignmentRead = 0;
	_previousAlignmentPos = -1;
};

void TBamFile::fill(TAlignment & alignment){
	//make sure container is empty
	alignment.clear();

	//add basic info
	//alignment.bamAlignment = bamAlignment; //DO WE NEED IT??

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

