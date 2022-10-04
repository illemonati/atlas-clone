/*
 * BamData.cpp
 *
 *  Created on: Jun 8, 2020
 *      Author: phaentu
 */

#include "TCigar.h"
#include <algorithm>
#include "stringFunctions.h"

namespace BAM {

//----------------------------------------------------------
// TCigar
// A class to store, access and manipulate CIGAR operators
//----------------------------------------------------------
TCigar::TCigar(TCigar cigar, uint16_t overlapLength, bool isFirst) {
	uint16_t overlap = 0;
	uint16_t nonOverlapLength = cigar.lengthRead() - overlapLength;
	if (isFirst) {
		auto iterator = cigar.begin();
		while (lengthRead() < nonOverlapLength) {
			std::cout << lengthRead() << " " << nonOverlapLength << std::endl;
			if (lengthRead()+iterator->length > nonOverlapLength) {
				if (iterator->type != 'S'){
					overlap = (lengthRead()+iterator->length) - nonOverlapLength;
					std::cout << "S:" << overlap << std::endl;
					add(iterator->type,iterator->length - overlap);
					iterator++;
				} else {
					overlap += iterator->length;
					nonOverlapLength = lengthRead();
					iterator++;
				}
			} else {
				add(iterator->type,iterator->length);
				iterator++;
			}
		}
		while (iterator != cigar.end()) {
			if (iterator->type == 'M' || iterator->type == 'I' || iterator->type == '=' || iterator->type == 'X' || iterator->type == 'S')
				overlap+=iterator->length;
			iterator++;
		}
		add('S',overlap);
	} else {
		auto iterator = --cigar.end();
		while (lengthRead() < nonOverlapLength) {
					if (lengthRead()+iterator->length > nonOverlapLength) {
						if (iterator->type != 'S'){
							overlap = (lengthRead()+iterator->length) - nonOverlapLength;
							add(iterator->type,iterator->length - overlap);
							iterator--;
						} else {
							overlap += iterator->length;
							nonOverlapLength = lengthRead();
							iterator--;
						}
					} else {
						add(iterator->type,iterator->length);
						iterator--;
					}
				}

		while (iterator != --cigar.begin()) {
			if (iterator->type == 'M' || iterator->type == 'I' || iterator->type == '=' || iterator->type == 'X' || iterator->type == 'S') {
				std::cout << iterator->type << iterator->length << std::endl;
				overlap+=iterator->length;
			}
			iterator--;
		}
		std::cout << overlap << std::endl;
		add('S',overlap);
		_flipCigar();
	}
}

void TCigar::_flipCigar() {
	std::reverse(_cigar.begin(), _cigar.end());
	uint32_t temp = lengthSoftClippedLeft();
	_lengthSoftClippedLeft = lengthSoftClippedRight();
	_lengthSoftClippedRight = temp;
}

void TCigar::clear() {
	_cigar.clear();
	_lengthAligned           = 0;
	_lengthInserted          = 0;
	_lengthDeleted           = 0;
	_lengthSkipped           = 0;
	_lengthSoftClippedLeft   = 0;
	_lengthSoftClippedRight  = 0;
};

void TCigar::add(char Type, uint32_t Length) {
	if (_lengthSoftClippedRight) { throw "Cigar string contains entries past soft clipping on right!"; }
	if (Type == 'M' || Type == '=' || Type == 'X') {
		_lengthAligned += Length;
	} else if (Type == 'I') {
		_lengthInserted += Length;
	} else if (Type == 'D') {
		_lengthDeleted += Length;
	} else if (Type == 'S') {
		if (!_lengthSoftClippedLeft) {
			_lengthSoftClippedLeft = Length;
		} else {
			_lengthSoftClippedRight = Length;
		}
	} else if (Type == 'N') {
		_lengthSkipped += Length;
	} else if (Type != 'H' && Type != 'P') {
		throw std::string("Unknown CIGAR operation '") + Type + "'!";
	}

	// add to vector
	_cigar.emplace_back(Type, Length);
};

void TCigar::removeSoftClips() {
	_cigar.erase(std::remove_if(_cigar.begin(), _cigar.end(), [](const auto &c) { return c.type == 'S'; }));

	// update length
	_lengthSoftClippedLeft  = 0;
	_lengthSoftClippedRight = 0;
};

std::string TCigar::compileString() const{
	std::string s;
	for(auto& it : _cigar){
		s += coretools::str::toString(it.length) + it.type;
	}
	return s;
}

TCigar::operator std::string() const {
	return compileString();
}

}; // namespace BAM
