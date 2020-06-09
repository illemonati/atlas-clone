/*
 * TFastaBuffer.h
 *
 *  Created on: Feb 10, 2020
 *      Author: wegmannd
 */

#ifndef TFASTABUFFER_H_
#define TFASTABUFFER_H_

#include "../bamtools/utils/bamtools_fasta.h"
#include "TGenotypeMap.h"

namespace BAM{

//-----------------------------------------------------
//TFastaBuffer
//-----------------------------------------------------
//a buffer class to speed up adding the reference sequence to each read
//This class makes use of the fact that bam files are sorted, hence the buffer can always start at the current position

class TFastaBuffer{
private:
	BamTools::Fasta _reference;
	uint32_t _bufferSize;
	std::string _referenceSequence;
	BAM::TGenomePosition _curPosition;
	int _curChr;
	long _curStart, _curEnd;
	bool _hasReference;

	void moveTo(const int & chr, const int32_t & pos);

public:
	TFastaBuffer(){
		_bufferSize = 100000;
		_referenceSequence = "";
		_curStart = -1;
		_curChr = -1;
		_curEnd = -1;
		_hasReference = false;
	};

	TFastaBuffer(std::string fastaFile){
		initialize(fastaFile);
	};

	TFastaBuffer(std::string fastaFile, const uint32_t BufferSize){
		initialize(fastaFile, BufferSize);
	};

	~TFastaBuffer(){};

	void initialize(std::string fastaFile, const uint32_t BufferSize=1000000);
	bool hasReference(){ return _hasReference; };
	void fill(const uint16_t & chr, const uint32_t start, const uint32_t end, std::string & ref);
	char refAt(const uint16_t & chr, const int32_t & position);
	int getCurChr(){return _curChr;};
	long getCurStart(){return _curStart;};
	long getCurEnd(){return _curEnd;};
};

}; //end namespace

#endif /* TFASTABUFFER_H_ */
