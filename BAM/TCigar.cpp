/*
 * BamData.cpp
 *
 *  Created on: Jun 8, 2020
 *      Author: phaentu
 */

#include "TCigar.h"
#include "coretools/devtools.h"
#include <algorithm>
#include "coretools/Strings/stringFunctions.h"

namespace BAM {

//----------------------------------------------------------
// TCigar
// A class to store, access and manipulate CIGAR operators
//----------------------------------------------------------
TCigar::TCigar(TCigar cigar, uint16_t overlapLength, bool isForwardStrand, size_t &mappedBasesClipped) {
	uint16_t overlap = 0;
	//if the overlap is larger than the number of aligned bases in the read (i.e. the other read eclipses this alignment), all bases are softclipped
	if (overlapLength >= cigar.lengthMapped()) {
		add('S', cigar.lengthRead());
	} else {
		//how many aligned bases before the overlap begins
		uint16_t nonOverlapLength = cigar.lengthMapped() - overlapLength;
		//if the read starts before its mate, go beginning->end, otherwise go end->beginning
		if (isForwardStrand) {
			std::vector<CigarOperator>::const_iterator iterator = cigar.begin();
			//copy cigar string until either the start of the overlap or the end of the read is reached
			int nextLengthMapped = 0;
			while (lengthMapped() < nonOverlapLength && iterator!=cigar.end()) {
				//if next segment of cigar string would exceed nonOverlapLength, split it up into M/=/X and S
				nextLengthMapped = lengthMapped()+iterator->length;
				if (nextLengthMapped >= nonOverlapLength && (iterator->type == 'D' || iterator->type == 'N')){
					nonOverlapLength = lengthMapped();
				} else if (nextLengthMapped > nonOverlapLength && (iterator->type =='M' || iterator->type == '=' || iterator->type == 'X')) {
						overlap = nextLengthMapped - nonOverlapLength;
						add(iterator->type,iterator->length - overlap);
				} else {
					add(iterator->type,iterator->length);
				}
				iterator++;
			}
			//sum up the length of the remaining overlap
			while (iterator != cigar.end()) {
				if (iterator->type == 'M' || iterator->type == 'I' || iterator->type == '=' || iterator->type == 'X' || iterator->type == 'S')
					overlap+=iterator->length;
				iterator++;
			}

			add('S',overlap);
		} else {
			//same as the part above, just use rbegin instead of begin to construct the cigar-string from right to left
			std::vector<CigarOperator>::const_reverse_iterator iterator = cigar.rbegin();
			int nextLengthMapped = 0;
			while (lengthMapped() < nonOverlapLength && iterator != cigar.rend()) {
				nextLengthMapped = lengthMapped()+iterator->length;
				if (nextLengthMapped >= nonOverlapLength && (iterator->type == 'D' || iterator->type == 'N')){
					mappedBasesClipped += iterator->length;
					nonOverlapLength = lengthMapped();
				} else if (nextLengthMapped > nonOverlapLength && (iterator->type =='M' || iterator->type == '=' || iterator->type == 'X')) {
						overlap = nextLengthMapped - nonOverlapLength;
						mappedBasesClipped += overlap;
						add(iterator->type,iterator->length - overlap);
				} else {
					add(iterator->type,iterator->length);
				}
				iterator++;
			}
			while (iterator != cigar.rend()) {
				if (iterator->type == 'M' ||  iterator->type == '=' || iterator->type == 'X') {
					mappedBasesClipped += iterator->length;
					overlap+=iterator->length;
				} else if (iterator->type == 'I' || iterator->type == 'S') {
					overlap+=iterator->length;
				} else if (iterator->type == 'D' || iterator->type == 'N') {
					mappedBasesClipped += iterator->length;
				}
				iterator++;
			}
			_lengthSoftClippedLeft++;
			add('S',overlap);
			_lengthSoftClippedLeft--;
			//cigar string needs to be flipped because it was created in reverse
			_flipCigar();
		}
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
	if (_lengthSoftClippedRight) { UERROR("Cigar string contains entries past soft clipping on right!"); }
	if (Type == 'M' || Type == '=' || Type == 'X') {
		_lengthAligned += Length;
	} else if (Type == 'I') {
		_lengthInserted += Length;
	} else if (Type == 'D') {
		_lengthDeleted += Length;
	} else if (Type == 'S') {
		if (_cigar.empty()) {
			_lengthSoftClippedLeft = Length;
		} else {
			_lengthSoftClippedRight = Length;
		}
	} else if (Type == 'N') {
		_lengthSkipped += Length;
	} else if (Type != 'H' && Type != 'P') {
		UERROR("Unknown CIGAR operation '", Type, "'!");
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

void TCigar::removeSoftClips(size_t maxNumberOfSoftClippedBases) {
	if(maxNumberOfSoftClippedBases == 0)
		_cigar.erase(std::remove_if(_cigar.begin(), _cigar.end(), [](const auto &c) { return c.type == 'S'; }));
	else {
		auto cIter = _cigar.begin(); 
		bool pastLeft = false;
		while (cIter != _cigar.end()) {
			if (cIter->type == 'S' && cIter->length > maxNumberOfSoftClippedBases){
				cIter->length = maxNumberOfSoftClippedBases;
				// if softclips are removed from the cigar string -> adjust internal softclip counts on the appropriate side
				pastLeft ? _lengthSoftClippedRight = maxNumberOfSoftClippedBases : _lengthSoftClippedLeft = maxNumberOfSoftClippedBases;
			}
			// as soon as mapped bases are reached, the iterator is past the softclips on the left
			if (cIter->type == 'M' || cIter->type == '=' || cIter->type == 'X')
				pastLeft = true;
			cIter++;
		}
	}
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
