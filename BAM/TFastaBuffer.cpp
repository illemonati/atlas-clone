/*
 * TFastaBuffer.cpp
 *
 *  Created on: Apr 2, 2020
 *      Author: vivian
 */

#include "TFastaBuffer.h"
#include <stdlib.h>
#include <algorithm>
#include <string>
#include "genometools/GenotypeTypes.h"

namespace BAM{

using genometools::Base;

void TFastaBuffer::_moveTo(const genometools::TGenomePosition Position) const{
	if(!_hasReference){
		throw "Can not move reference: no FASTA file provided!";
	}

	_coordinates.move(Position, Position + _bufferSize);

	int fastaLength;
	if (!_reference.GetLength(_coordinates.refID(), fastaLength))
		throw "Problem reading Fasta length " + std::to_string(_coordinates.refID()) + " from index file!";

	const int start          = _coordinates.fromOnChr();
	const int end            = _coordinates.toOnChr();
	const std::string before = start < 0 ? std::string(std::abs(start), 'N') : "";
	const std::string after  = end > fastaLength ? std::string(end - fastaLength, 'N') : "";

	if(!_reference.GetSequence(_coordinates.refID(), std::max(0, start), std::min(end, fastaLength), _referenceSequence))
		throw "Problem reading " + std::to_string(_coordinates.refID()) + ":" + std::to_string(_coordinates.fromOnChr()) + "-" + std::to_string(_coordinates.toOnChr()) + " from fasta file! Are you using the correct fasta file?";
	_referenceSequence = before + _referenceSequence + after;
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

void TFastaBuffer::_fill(std::vector<Base> & VecToFill, size_t Length, size_t Offset) const {
	VecToFill.resize(Length);
	for(size_t p = 0; p < Length; ++p){
		VecToFill[p] = genometools::char2base(_referenceSequence[Offset + p]);
	}
};

void TFastaBuffer::fill(const genometools::TGenomeWindow & Window, std::vector<Base> & VecToFill) const {
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

void TFastaBuffer::fill(const genometools::TGenomePosition & Position, uint32_t Length, std::vector<Base> & VecToFill) const {
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

void TFastaBuffer::fill(const genometools::TGenomePosition & Start, const genometools::TGenomePosition & End, std::vector<Base> & VecToFill) const {
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

char TFastaBuffer::refCharAt(const genometools::TGenomePosition Position) const {
	if(Position  < _coordinates || Position > _coordinates){
		_moveTo(Position);
	}
	return _referenceSequence[Position - _coordinates.from()];
};

Base TFastaBuffer::refAt(const genometools::TGenomePosition Position) const{
	return genometools::char2base(refCharAt(Position));
};

}; //end namesapce
