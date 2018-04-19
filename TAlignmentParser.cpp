/*
 * TAlignmentParser.cpp
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#include "TAlignmentParser.h"

//-----------------------------------------------------
//TFastaBuffer
//-----------------------------------------------------
TFastaBuffer::TFastaBuffer(BamTools::Fasta* Reference){
	bufferSize = 10000000;
	reference = Reference;
	referenceSequence = "";
	curStart = -1;
	curChr = -1;
	curEnd = -1;
}


void TFastaBuffer::moveTo(const int & chr, const int32_t & pos){
	curChr = chr;
	curStart = pos;
	curEnd = pos + bufferSize;
	if(!reference->GetSequence(chr, curStart, curEnd, referenceSequence))
		throw "Problem reading " + toString(chr) + ":" + toString(curStart) + "-" + toString(curEnd) + " from fasta file!";
}

void TFastaBuffer::fill(const int & chr, const int32_t & start, const int32_t end, std::string & ref){
	//move buffer, if necessary
	if(chr != curChr || end > curEnd || start < curStart)
		moveTo(chr, start);

	//now copy to string
	ref.assign(referenceSequence, start - curStart, end - start + 1);
}

//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
TAlignmentParser::TAlignmentParser(){
	maxSize = 0;
	readGroupTable = NULL;
	logfile = NULL;
	_keepDuplicates = false;
	initialized = false;

	//filters
	applyQualityFilter = false;
	minQual = 33;
	maxQual = 126;
	minQualForPrinting = 33;
	maxQualForPrinting = 126;
	trimReads = false;
	trimmingLength3Prime = 0;
	trimmingLength5Prime = 0;

	//reference
	hasReference = false;
	fastaBuffer = NULL;
};

TAlignmentParser::TAlignmentParser(TReadGroups* ReadGroupTable, unsigned int MaxSize, TLog* Logfile){
	TAlignmentParser();
	init(ReadGroupTable, MaxSize, Logfile);
};

void TAlignmentParser::init(TReadGroups* ReadGroupTable, unsigned int MaxSize, TLog* Logfile){
	clear();

	maxSize = MaxSize;
	readGroupTable = ReadGroupTable;
	logfile = Logfile;

	//data
	base = new Base[maxSize];
	baseAsChar = new char[maxSize];
	context = new BaseContext[maxSize];
	qualityOriginal = new int[maxSize];
	qualityRecalibrated = new int[maxSize];
	errorRates = new double[maxSize];
	aligned = new bool[maxSize];
	alignedPos = new int[maxSize];
	distFrom3Prime = new int[maxSize];
	distFrom5Prime = new int[maxSize];
	pmdCT = new double[maxSize];
	pmdGA = new double[maxSize];

	//soft clipped data
	softClippedLength = new int[2];
	softClippedBase = new char*[2];
	softClippedQuality = new char*[2];
	softClippedBase[0] = new char[maxSize];
	softClippedBase[1] = new char[maxSize];
	softClippedQuality[0] = new char[maxSize];
	softClippedQuality[1] = new char[maxSize];

	initialized = true;
};

void TAlignmentParser::clear(){
	if(initialized){
		delete[] base;
		delete[] baseAsChar;
		delete[] context;
		delete[] qualityOriginal;
		delete[] qualityRecalibrated;
		delete[] errorRates;
		delete[] aligned;
		delete[] alignedPos;
		delete[] distFrom3Prime;
		delete[] distFrom5Prime;
		delete[] pmdCT;
		delete[] pmdGA;

		delete[] softClippedLength;
		delete[] softClippedBase[0];
		delete[] softClippedBase[1];
		delete[] softClippedBase;
		delete[] softClippedQuality[0];
		delete[] softClippedQuality[1];
		delete[] softClippedQuality;

		initialized = false;
	}
};

void TAlignmentParser::setQualityFilters(int MinQual, int MaxQual){
	applyQualityFilter = true;
	minQual = MinQual;
	maxQual = MaxQual;
};

void TAlignmentParser::setQualityRangeForPrinting(int minQual, int maxQual){
	minQualForPrinting = minQual;
	maxQualForPrinting = maxQual;
};

void TAlignmentParser::setReadTrimming(int trim3Prime, int trim5Prime){
	trimmingLength3Prime = trim3Prime;
	trimmingLength5Prime = trim5Prime;
	trimReads = true;
};

void TAlignmentParser::addReference(BamTools::Fasta* reference){
	hasReference = true;
	fastaBuffer = new TFastaBuffer(reference);
};




//------------------------------
//public functions
//------------------------------
bool TAlignmentParser::readAlignment(BamTools::BamReader & bamReader, TAlignment & alignment){
	if(!bamReader.GetNextAlignment(bamAlignment))
		return false;

	//check read length
	if(bamAlignment.AlignedBases.size() > maxSize)
		throw "Alignment '" +  bamAlignment.Name + "' is longer than the max read length! Please change max read length to parse this data.";

	//fill alignment
	std::string readGroup;
	bamAlignment.GetTag("RG", readGroup);
	int readGroupId = readGroupTable->find(readGroup);
	alignment.fill(bamAlignment, readGroupId);

	//check if insert size is shorter than read, this means we are reading the adaptor sequence
	bool filtersPassed = true;
	if(bamAlignment.IsPaired() && abs(bamAlignment.InsertSize) <= bamAlignment.AlignedBases.length()){
		logfile->warning("The following alignment is longer than its insert size: " + bamAlignment.Name);
		filtersPassed = false;
	} else {
		//apply filters: read group in use and basic QC
		filtersPassed = readGroupTable->readGroupInUse(readGroupId)
						&& bamAlignment.IsMapped() && !bamAlignment.IsFailedQC()
						&& bamAlignment.IsPrimaryAlignment()
						&& (_keepDuplicates || !bamAlignment.IsDuplicate());
	}

	//set filter
	alignment.setFiltersPassed(filtersPassed);

	return true;
}







