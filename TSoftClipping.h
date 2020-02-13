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

	void empty();
	void assessSoftClipping(BamTools::BamAlignment bamAlignment);
};

//--------------------------------------------------------
// TSoftClippingStatsFile
//--------------------------------------------------------
class TSoftClippingStatsFile{
private:
	TOutputFileZipped out;
	bool printSoftClippedSequences;


public:
	TSoftClippingStatsFile(const std::string outputName, bool PrintSoftClippedSequences);
	void write(const std::string & readName, const int & position, TSoftClippingData & data);
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
	TSoftClippingMatrixStorage();
	TSoftClippingMatrixStorage(int MaxReadLength);
	~TSoftClippingMatrixStorage(){ clear(); };
	void allocate(int MaxReadLength);
	void clear();
	void empty();
	void add(int readLength, int softClippedLength);
	long size();
	void write(std::string filename);
};

//--------------------------------------------------------
// TSoftClippingMatrix
// Counts read length and soft clipping length
//--------------------------------------------------------
class TSoftClippingMatrix{
private:
	TSoftClippingMatrixStorage _left, _right, _total;

public:
	TSoftClippingMatrix(int MaxReadLength);
	~TSoftClippingMatrix(){};

	void clear();
	void add(const TSoftClippingData & data);
	void write(const std::string & outputName);
};


#endif /* TSOFTCLIPPING_H_ */
