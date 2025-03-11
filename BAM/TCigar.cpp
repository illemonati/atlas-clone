/*
 * BamData.cpp
 *
 *  Created on: Jun 8, 2020
 *      Author: phaentu
 */

#include "TCigar.h"
#include "coretools/Main/TError.h"
#include <algorithm>

namespace BAM {

//----------------------------------------------------------
// TCigar
// A class to store, access and manipulate CIGAR operators
//----------------------------------------------------------

void TCigar::add(char Type, size_t Length) {
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
}

void TCigar::_compileLengths() {
	_lengthInserted = 0;
	_lengthDeleted  = 0;
	_lengthAligned  = 0;
	_lengthSkipped  = 0;

	_lengthSoftClippedLeft  = 0;
	_lengthSoftClippedRight = 0;
	for (const auto &c : _cigar) {
		switch (c.type) {
		case 'I': _lengthInserted += c.length; break;
		case 'D':_lengthDeleted += c.length; break;
		case 'N': _lengthSkipped += c.length; break;
		case 'M':
		case '=':
		case 'X': _lengthAligned += c.length; break;
		case 'S': break;
		default: DEVERROR("Error parsing cigar '", compileString(), "'.");
		}
	}
	if (_cigar.front().type == 'S') _lengthSoftClippedLeft = _cigar.front().length;
	if (_cigar.size() > 1 && _cigar.back().type == 'S') _lengthSoftClippedLeft = _cigar.back().length;
}

void TCigar::removeSoftClips() {
	if (_cigar.front().type == 'S') _cigar.erase(_cigar.begin());
	if (_cigar.back().type == 'S') _cigar.pop_back();

	// update length
	_lengthSoftClippedLeft  = 0;
	_lengthSoftClippedRight = 0;
}

size_t TCigar::removeMappedLeft(size_t Length) {
	std::reverse(_cigar.begin(), _cigar.end());
	removeMappedRight(Length);
	std::reverse(_cigar.begin(), _cigar.end());

	std::swap(_lengthSoftClippedLeft, _lengthSoftClippedRight);
	return Length;
}

void TCigar::removeMappedRight(size_t Length) {
	if (Length == 0) return;

	if (Length > lengthMapped())
		DEVERROR("Cannot add ", Length, " Softclips to cigar '", compileString(), "'.");

	CigarOperator softClipR('S', 0);
	if (_cigar.back().type == 'S') {
		softClipR.length += _cigar.back().length;
		_cigar.pop_back();
	}

	while (Length > 0 && !_cigar.empty()) {
		switch (_cigar.back().type) {
		case 'I':
			softClipR.length += _cigar.back().length;
			_cigar.pop_back();
			break;
		case 'D':
		case 'N':
			if (Length >= _cigar.back().length) {
				Length -= _cigar.back().length;
				_cigar.pop_back();
			} else {
				_cigar.back().length -= Length;
				Length = 0;
			}
			break;

		case 'M':
		case '=':
		case 'X':
			if (Length >= _cigar.back().length) {
				Length -= _cigar.back().length;
				softClipR.length += _cigar.back().length;
				_cigar.pop_back();
			} else {
				_cigar.back().length -= Length;
				softClipR.length += Length;
				Length = 0;
			}
			break;
		default: DEVERROR("Error parsing cigar '", compileString(), "'.");
		}
	}
	if (Length > 0) DEVERROR("Error parsing cigar '", compileString(), "'.");
	_cigar.push_back(softClipR);

	if (_cigar.size() == 2 && _cigar.front().type == 'S' && _cigar.back().type == 'S') {
		// merge softClips
		_cigar.front().length += _cigar.back().length;
		_cigar.pop_back();
	}
	_compileLengths();
}

void TCigar::trimSoftClips(size_t maxNumberOfSoftClippedBases) {
	if(maxNumberOfSoftClippedBases == 0) removeSoftClips();
	else {
		if (_cigar.front().type == 'S' && _cigar.front().length > maxNumberOfSoftClippedBases) {
			_cigar.front().length = maxNumberOfSoftClippedBases;
			_lengthSoftClippedLeft = maxNumberOfSoftClippedBases;
		}
		if (_cigar.back().type == 'S' && _cigar.back().length > maxNumberOfSoftClippedBases) {
			_cigar.back().length = maxNumberOfSoftClippedBases;
			_lengthSoftClippedRight = maxNumberOfSoftClippedBases;
		}
	}
}

std::string TCigar::compileString() const{
	std::string s;
	for(const auto& c : _cigar){
		s += coretools::str::toString(c.length) + c.type;
	}
	return s;
}

TCigar::operator std::string() const {
	return compileString();
}

void TCigar::clear() {
	_cigar.clear();
	_lengthAligned          = 0;
	_lengthInserted         = 0;
	_lengthDeleted          = 0;
	_lengthSkipped          = 0;
	_lengthSoftClippedLeft  = 0;
	_lengthSoftClippedRight = 0;
}

} // namespace BAM
