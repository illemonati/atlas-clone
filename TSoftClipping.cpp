/*
 * TSoftClipping.cpp
 *
 *  Created on: Dec 12, 2019
 *      Author: schulzi
 */

#include "TSoftClipping.h"

//------------------------------------------
// TSoftClippingData
//------------------------------------------
void TSoftClippingData::empty(){
	hasSoftClipping = false;
	alignmentLength = 0;
	softClippingLength_left = 0;
	softClippingLength_right = 0;
	middleLength = 0;
	softClippingBases_left.clear();
	softClippingBases_right.clear();
	middleBases.clear();
	middleQualities.clear();
};

void TSoftClippingData::assessSoftClipping(BamTools::BamAlignment bamAlignment){
	//count S, not S, S pattern from cigar string
	empty();
	bool reachedMiddle = false;

	std::vector<BamTools::CigarOp>::const_iterator cigarIter = bamAlignment.CigarData.begin();
	std::vector<BamTools::CigarOp>::const_iterator cigarEnd  = bamAlignment.CigarData.end();

	//position in read, i is position in cigar type
	int p = 0;

	for(; cigarIter != cigarEnd; ++cigarIter ){
		if(cigarIter->Type == 'D')
			continue;
		if(cigarIter->Type == 'S'){
			if(reachedMiddle){
				softClippingLength_right += cigarIter->Length;
				softClippingBases_right.append(bamAlignment.QueryBases, p, cigarIter->Length);
			} else {
				softClippingLength_left += cigarIter->Length;
				softClippingBases_left.append(bamAlignment.QueryBases, p, cigarIter->Length);

			}
		} else {
			reachedMiddle = true;
			middleLength += cigarIter->Length;
			middleBases.append(bamAlignment.QueryBases, p, cigarIter->Length);
			middleQualities.append(bamAlignment.Qualities, p, cigarIter->Length);
		}
		p += cigarIter->Length;
	}

	//set other things
	alignmentLength = bamAlignment.QueryBases.length();
	softClippingLength = softClippingLength_left + softClippingLength_right;
	if(softClippingLength_left + softClippingLength_right > 0)
		hasSoftClipping = true;


	//std::cout << "alignment '" << bamAlignment.Name << "':";
	//std::cout << "Softclipping (" << hasSoftClipping << ") length = " << softClippingLength << " = " << softClippingLength_left << " + " << softClippingLength_right << std::endl;
};

//------------------------------------------
// TSoftClippingStatsFile
//------------------------------------------
TSoftClippingStatsFile::TSoftClippingStatsFile(const std::string outputName, bool PrintSoftClippedSequences){
	printSoftClippedSequences = PrintSoftClippedSequences;

	//open output file
	out.open(outputName + "_softClippingStats.txt.gz");

	//write header
	std::vector<std::string> header = {"Read", "position", "nClippedLeft", "nNotClipped", "nClippedRight"};
	if(printSoftClippedSequences)
		header.push_back({"sNotClipped", "sClippedRight"});
	out.writeHeader(header);
};

void TSoftClippingStatsFile::write(const std::string & readName, const int & position, TSoftClippingData & data){
	out << readName << position << data.softClippingLength_left << data.middleLength << data.softClippingLength_right;
	if(printSoftClippedSequences)
		out << data.softClippingBases_left << data.softClippingBases_right;
	out << std::endl;
};


//--------------------------------------------------------
// TSoftClippingMatrixStorage
//--------------------------------------------------------
TSoftClippingMatrixStorage::TSoftClippingMatrixStorage(){
	_allocated = false;
	_maxReadLength = 0;
	_counts = NULL;
};

TSoftClippingMatrixStorage::TSoftClippingMatrixStorage(int MaxReadLength){
	_allocated = false;
	allocate(MaxReadLength);
};

void TSoftClippingMatrixStorage::allocate(int MaxReadLength){
	_maxReadLength = MaxReadLength;
	if(_allocated) clear();
	_counts = new int*[_maxReadLength];
	for(int i=0; i<_maxReadLength; ++i){
		_counts[i] = new int[i+1];
	}
	_allocated = true;

	empty();
};

void TSoftClippingMatrixStorage::clear(){
	for(int i=0; i<_maxReadLength; ++i){
		delete[] _counts[i];
	}
	delete[] _counts;
	_allocated = false;
};

void TSoftClippingMatrixStorage::empty(){
	for(int i=0; i<_maxReadLength; ++i){
		for(int j=0; j<(i+1); ++j){
			_counts[i][j]=0;
		}
	}
};

void TSoftClippingMatrixStorage::add(int readLength, int softClippedLength){
	++_counts[readLength][softClippedLength];
};

long TSoftClippingMatrixStorage::size(){
	long sum = 0;
	for(int i=0; i<_maxReadLength; ++i){
		for(int j=0; j<(i+1); ++j){
			sum += _counts[i][j];
		}
	}
	return sum;
};

void TSoftClippingMatrixStorage::write(std::string filename){
	//find longest readlength with counts
	int maxLengthForWriting = 0;
	for(int i=0; i<_maxReadLength; ++i){
		for(int j=0; j<(i+1); ++j){
			if(_counts[i][j] > 0 && i > maxLengthForWriting)
				maxLengthForWriting = i;
		}
	}
	++maxLengthForWriting;

	//open table file
	TOutputFilePlain out(filename);

	//write header
	std::vector<std::string> header = {"readLength"};
	for(int i=0; i <= maxLengthForWriting; ++i){
		header.push_back("softClippedLength_" + toString(i));
	}
	out.writeHeader(header);

	//write counts
	for(int i=1; i < maxLengthForWriting; ++i){
		out << i;

		for(int j=0; j<(i+1); ++j){
			out << _counts[i][j];
		}

		for(int j=i+1; j <= maxLengthForWriting; ++j){
			out << "0";
		}

		out << std::endl;
	}
};

//--------------------------------------------------------
// TSoftClippingMatrix
// Counts read length and soft clipping length
//--------------------------------------------------------
TSoftClippingMatrix::TSoftClippingMatrix(int MaxReadLength){
	//allocate matrices
	_left.allocate(MaxReadLength);
	_right.allocate(MaxReadLength);
	_total.allocate(MaxReadLength);
};

void TSoftClippingMatrix::clear(){
	_left.clear();
	_right.clear();
	_total.clear();
};


void TSoftClippingMatrix::add(const TSoftClippingData & data){
	_left.add(data.alignmentLength, data.softClippingLength_left);
	_right.add(data.alignmentLength, data.softClippingLength_right);
	_total.add(data.alignmentLength, data.softClippingLength);
};

void TSoftClippingMatrix::write(const std::string & outputName){
	_left.write(outputName + "_softClippingMatrixLeft.txt");
	_right.write(outputName + "_softClippingMatrixRight.txt");
	_total.write(outputName + "_softClippingMatrixTotal.txt");
};


