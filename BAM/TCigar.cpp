/*
 * BamData.cpp
 *
 *  Created on: Jun 8, 2020
 *      Author: phaentu
 */

#include "TCigar.h"
#include <algorithm>
#include <string>

namespace BAM {

//-----------------------------------------------------
// TCigar
// A class to store, access and manipulate CIGAR operators
//-----------------------------------------------------
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
	if (_lengthSoftClippedRight) { throw "Cigar string contains entries past softclipping on right!"; }
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

}; // namespace BAM
