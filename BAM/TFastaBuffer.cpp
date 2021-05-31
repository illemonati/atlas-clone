/*
 * TFastaBuffer.cpp
 *
 *  Created on: Apr 2, 2020
 *      Author: vivian
 */


#include "TFastaBuffer.h"

namespace BAM{

void TFastaBuffer::_moveTo(const TGenomePosition Position) const{
	//NOTE: bamtools was modified to append N in case pos < 0 or pos+length > is beyond chromosome. This is the expected behavior in ATLAS and must be preserved!
	if(!_hasReference){
		throw "Can not move reference: no FASTA file provided!";
	}

	_coordinates.move(Position, Position + _bufferSize);

	if(!_reference.GetSequence(_coordinates.refID(), _coordinates.fromOnChr(), _coordinates.toOnChr(), _referenceSequence))
		throw "Problem reading " + toString(_coordinates.refID()) + ":" + toString(_coordinates.fromOnChr()) + "-" + toString(_coordinates.toOnChr()) + " from fasta file! Are you using the correct fasta file?";
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
	_hasReference = true;
};

void TFastaBuffer::_fill(std::vector<Base> & VecToFill, const size_t & Length, const size_t & Offset) const {
	VecToFill.resize(Length);
	for(size_t p = 0; p < Length; ++p){
		VecToFill[p] = _referenceSequence[Offset + p];
	}
};

void TFastaBuffer::fill(const TGenomeWindow & Window, std::vector<Base> & VecToFill) const {
	//move buffer, if necessary
	if(!_coordinates.contains(Window)){
		//adjust buffer size?
		if(Window.size() > _bufferSize){
			_bufferSize = Window.size();
		}

		//move
		_moveTo(Window.from());
	}

	//fill
	_fill(VecToFill, Window.size(), Window - _coordinates);

};

void TFastaBuffer::fill(const TGenomePosition & Position, const uint32_t & Length, std::vector<Base> & VecToFill) const {
	//move buffer, if necessary
	if(Position  < _coordinates || Position + Length - 1 > _coordinates){
		if(Length > _bufferSize){
			_bufferSize = Length;
		}
		_moveTo(Position);
	}

	//fill
	_fill(VecToFill, Length, Position - _coordinates);
};

void TFastaBuffer::fill(const TGenomePosition & Start, const TGenomePosition & End, std::vector<Base> & VecToFill) const {
	fill(Start, End - Start + 1, VecToFill);
};

/*
void TFastaBuffer::fill(const TGenomePosition & Position, const uint32_t end, std::string & ref){
	//move buffer, if necessary
	if(Pos)
	if(chr != _curChr || end > _curEnd || start < _curStart){
		if(end - start + 1 > _bufferSize){
			_bufferSize = end - start + 1;
		}
		_moveTo(chr, start);
	}

	//now copy to string
	ref.assign(_referenceSequence, start - _curStart, end - start + 1);
};
*/

char TFastaBuffer::refCharAt(const TGenomePosition Position) const {
	if(Position  < _coordinates || Position > _coordinates){
		_moveTo(Position);
	}
	return _referenceSequence[Position - _coordinates.from()];
};

Base TFastaBuffer::refAt(const TGenomePosition Position) const{
	return Base(refCharAt(Position));
};

}; //end namesapce
