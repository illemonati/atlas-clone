/*
 * TFastaBuffer.cpp
 *
 *  Created on: Apr 2, 2020
 *      Author: vivian
 */


#include "TFastaBuffer.h"

void TFastaBuffer::moveTo(const int & chr, const int32_t & pos){
	curChr = chr;
	curStart = pos;
	curEnd = pos + bufferSize;
	if(!reference->GetSequence(chr, curStart, curEnd, referenceSequence))
		throw "Problem reading " + toString(chr) + ":" + toString(curStart) + "-" + toString(curEnd) + " from fasta file!";
};

void TFastaBuffer::initialize(BamTools::Fasta* Reference, const uint32_t BufferSize){
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

void TFastaBuffer::fill(const int & chr, const uint32_t start, const uint32_t end, std::string & ref){
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

char TFastaBuffer::refAt(const int & chr, const int32_t & position){
	if(chr != curChr || position > curEnd || position < curStart){
		moveTo(chr, position);
	}
	return referenceSequence[position - curStart];
};
