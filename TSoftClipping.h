/*
 * TSoftClipping.h
 *
 *  Created on: Dec 9, 2019
 *      Author: wegmannd
 */

#ifndef TSOFTCLIPPING_H_
#define TSOFTCLIPPING_H_

#include "bamtools/api/BamAlignment.h"
#include "TFile.h"

//--------------------------------------------------------
// TSoftClippingData
//--------------------------------------------------------
class TSoftClippingData{
public:
	bool hasSoftClipping;
	int alignmentLength;
	int softClippingLength, softClippingLength_left, softClippingLength_right;
	int middleLength; //Everything in between softclipped bases
	std::string softClippingBases_left, softClippingBases_right;
	std::string middleBases, middleQualities;

	TSoftClippingData(){
		empty();
	};

	void empty(){
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

	void assessSoftClipping(BamTools::BamAlignment bamAlignment){
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
	};
};

//--------------------------------------------------------
// TSoftClippingStatsFile
//--------------------------------------------------------
class TSoftClippingStatsFile{
private:
	TOutputFileZipped out;
	bool printSoftClippedSequences;


public:
	TSoftClippingStatsFile(const std::string outputName, bool PrintSoftClippedSequences){
		printSoftClippedSequences = PrintSoftClippedSequences;

		//open output file
		out.open(outputName + "_softClippingStats.txt.gz");

		//write header
		std::vector<std::string> header = {"Read", "position", "nClippedLeft", "nNotClipped", "nClippedRight"};
		if(printSoftClippedSequences)
			header.push_back({"sNotClipped", "sClippedRight"});
		out.writeHeader(header);
	};

	void write(const std::string & readName, const int & position, TSoftClippingData & data){
		out << readName << position << data.softClippingLength_left, data.middleLength, data.softClippingLength_right;
		if(printSoftClippedSequences)
			out << data.softClippingBases_left << data.softClippingBases_right;
		out << std::endl;
	};
};

//--------------------------------------------------------
// TSoftClippingMatrixStorage
//--------------------------------------------------------
class TSoftClippingMatrixStorage{
private:
	int** _counts; // [readLength][softClippedLength]
	bool _allocated;
	int _maxReadLength;

public:
	TSoftClippingMatrixStorage(){
		_allocated = false;
		_maxReadLength = 0;
		_counts = NULL;
	};

	TSoftClippingMatrixStorage(int MaxReadLength){
		_allocated = false;
		allocate(MaxReadLength);
	};

	void allocate(int MaxReadLength){
		_maxReadLength = MaxReadLength;
		if(_allocated) clear();
		_counts = new int*[_maxReadLength];
		for(int i=0; i<_maxReadLength; ++i){
			_counts[i] = new int[i+1];
	    }
		_allocated = true;

		empty();
	};

	void clear(){
		for(int i=0; i<_maxReadLength; ++i){
			delete[] _counts[i];
		}
		delete[] _counts;
		_allocated = false;
	};

	void empty(){
		for(int i=0; i<_maxReadLength; ++i){
			for(int j=0; j<(i+1); ++j){
				_counts[i][j]=0;
			}
		}
	};

	void add(int readLength, int softClippedLength){
		++_counts[readLength][softClippedLength];
	};

	void write(std::string filename){
		//open table file
		TOutputFilePlain out(filename);

		//write header
		std::vector<std::string> header = {"readLength / softClippedLength"};
		for(int i=0; i <= _maxReadLength; ++i){
			header.push_back(toString(i));
		}
		out.writeHeader(header);

		//write counts
		for(int i=0; i < _maxReadLength; ++i){
			out << i+1;
			for(int j=0; j<(i+1); ++j){
				out << _counts[i][j];
			}
			for(int j=i+1; j <= _maxReadLength; ++j){
				out << "0";
			}
			out << std::endl;
		}
	};
};

//--------------------------------------------------------
// TSoftClippingMatrix
// Counts read length and soft clipping length
//--------------------------------------------------------
class TSoftClippingMatrix{
private:
	TSoftClippingMatrixStorage _left, _right, _total;

public:
	TSoftClippingMatrix(int MaxReadLength){
		//allocate matrices
		_left.allocate(MaxReadLength);
		_right.allocate(MaxReadLength);
		_total.allocate(MaxReadLength);
	};

	~TSoftClippingMatrix(){};

	void clear(){
	};

	void add(const TSoftClippingData & data){
		_left.add(data.alignmentLength, data.softClippingLength_left);
		_right.add(data.alignmentLength, data.softClippingLength_right);
		_total.add(data.alignmentLength, data.softClippingLength);
	};

	void write(const std::string & outputName){
		_left.write(outputName + "_softClippingMatrixLeft.txt");
		_right.write(outputName + "_softClippingMatrixRight.txt");
		_total.write(outputName + "_softClippingMatrixTotal.txt");
	};
};


#endif /* TSOFTCLIPPING_H_ */
