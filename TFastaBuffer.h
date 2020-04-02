/*
 * TFastaBuffer.h
 *
 *  Created on: Feb 10, 2020
 *      Author: wegmannd
 */

#ifndef TFASTABUFFER_H_
#define TFASTABUFFER_H_

#include "bamtools/utils/bamtools_fasta.h"
#include "TGenotypeMap.h"

//-----------------------------------------------------
//TFastaBuffer
//-----------------------------------------------------
//a buffer class to speed up adding the reference sequence to each read
//This class makes use of the fact that bam files are sorted, hence the buffer can always start at the current position

class TFastaBuffer{
private:
	BamTools::Fasta* reference;
	uint32_t bufferSize;
	std::string referenceSequence;
	int curChr;
	long curStart, curEnd;

	void moveTo(const int & chr, const int32_t & pos);

public:
	TFastaBuffer(){
		bufferSize = 100000;
		reference = nullptr;
		referenceSequence = "";
		curStart = -1;
		curChr = -1;
		curEnd = -1;
	};

	TFastaBuffer(BamTools::Fasta* Reference){
		initialize(Reference, 1000000);
	};

	TFastaBuffer(BamTools::Fasta* Reference, const uint32_t BufferSize){
		initialize(Reference, BufferSize);
	};

	~TFastaBuffer(){};

	void initialize(BamTools::Fasta* Reference, const uint32_t BufferSize);
	void fill(const int & chr, const uint32_t start, const uint32_t end, std::string & ref);
	char refAt(const int & chr, const int32_t & position);
	int getCurChr(){return curChr;};
	long getCurStart(){return curStart;};
	long getCurEnd(){return curEnd;};
};



#endif /* TFASTABUFFER_H_ */
