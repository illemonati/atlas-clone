/*
 * TFastaBuffer.cpp
 *
 *  Created on: Apr 2, 2020
 *      Author: vivian
 */


#include "TFastaBuffer.h"

namespace BAM{

void TFastaBuffer::moveTo(const int & chr, const int32_t & pos){
	//NOTE: bamtools was modified to append N in case pos < 0 or pos+length > is beyond chromosome. This is the expected behavior in ATLAS and must be preserved!
	if(!_hasReference){
		throw "Can not move reference: no FASTA file provided!";
	}
	_curChr = chr;
	_curStart = pos;
	_curEnd = pos + _bufferSize;
	if(!_reference.GetSequence(chr, _curStart, _curEnd, _referenceSequence))
		throw "Problem reading " + toString(chr) + ":" + toString(_curStart) + "-" + toString(_curEnd) + " from fasta file!";
};

void TFastaBuffer::initialize(std::string fastaFile, const uint32_t BufferSize){
	if(BufferSize < 1000)
		_bufferSize = 1000;
	else
		_bufferSize = BufferSize;

	//open fasta reference
	std::string fastaIndex = fastaFile + ".fai";
	if(!_reference.Open(fastaFile, fastaIndex))
		throw "Failed to open FASTA file '" + fastaFile + "'! Is index file present?";

	//initialize tmp variables
	_referenceSequence = "";
	_curStart = -1;
	_curChr = -1;
	_curEnd = -1;
};

void TFastaBuffer::fill(const uint16_t & chr, const uint32_t start, const uint32_t end, std::string & ref){
	//move buffer, if necessary
	if(chr != _curChr || end > _curEnd || start < _curStart){
		if(end - start + 1 > _bufferSize){
			_bufferSize = end - start + 1;
		}
		moveTo(chr, start);
	}

	//now copy to string
	ref.assign(_referenceSequence, start - _curStart, end - start + 1);
};

char TFastaBuffer::refAt(const uint16_t & chr, const int32_t & position){
	if(chr != _curChr || position > _curEnd || position < _curStart){
		moveTo(chr, position);
	}
	return _referenceSequence[position - _curStart];
};

}; //end namesapce
