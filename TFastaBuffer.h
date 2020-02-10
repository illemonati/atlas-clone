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

	void moveTo(const int & chr, const int32_t & pos){
		curChr = chr;
		curStart = pos;
		curEnd = pos + bufferSize;
		if(!reference->GetSequence(chr, curStart, curEnd, referenceSequence))
			throw "Problem reading " + toString(chr) + ":" + toString(curStart) + "-" + toString(curEnd) + " from fasta file!";
	};


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

	void initialize(BamTools::Fasta* Reference, const uint32_t BufferSize){
		if(BufferSize < 1000)
			bufferSize = 1000;
		else
			bufferSize = BufferSize;
		reference = Reference;
		referenceSequence = "";
		curStart = -1;
		curChr = -1;
		curEnd = -1;
	};

	void fill(const int & chr, const uint32_t start, const uint32_t end, std::string & ref){
		//move buffer, if necessary
		if(chr != curChr || end > curEnd || start < curStart){
			if(end - start + 1 > bufferSize){
				bufferSize = end - start + 1;
			}
			moveTo(chr, start);
		}

		//now copy to string
		ref.assign(referenceSequence, start - curStart, end - start + 1);
	};

	char refAt(const int & chr, const int32_t & position){
		if(chr != curChr || position > curEnd || position < curStart){
			moveTo(chr, position);
		}
		return referenceSequence[position - curStart];
	};

	int getCurChr(){return curChr;};
	long getCurStart(){return curStart;};
	long getCurEnd(){return curEnd;};
};



#endif /* TFASTABUFFER_H_ */
