/*
 * BamData.cpp
 *
 *  Created on: Jun 8, 2020
 *      Author: phaentu
 */

#include "TCigar.h"

namespace BAM{

//-----------------------------------------------------
//TCigar
//A class to store, access and manipulate CIGAR operators
//-----------------------------------------------------
void TCigar::clear(){
	_cigar.clear();
	_lengthAligned = 0;
	_lengthInserted = 0;
	_lengthDeleted = 0;
	_lengthSoftClippedLeft = 0;
	_lengthSoftClippedRight = 0;
	_addSoftClippedLeft = true;
};

void TCigar::add(const char & Type, const uint32_t & Length){
	if(Type == 'M' || Type == '=' || Type == 'X'){
		_lengthAligned += Length;
		_addSoftClippedLeft = false;
	} else if(Type == 'I'){
		_lengthInserted += Length;
		_addSoftClippedLeft = false;
	} else if(Type =='D'){
		_lengthDeleted += Length;
		_addSoftClippedLeft = false;
	} else if(Type == 'S'){
		if(_addSoftClippedLeft){
			_lengthSoftClippedLeft += Length;
		} else {
			_lengthSoftClippedRight += Length;
		}
	} else if(Type != 'N' && Type != 'H' && Type != 'P'){
		throw "Unknown CIGAR operation '" + Type + "'!";
	}

	//add to vector
	_cigar.emplace_back(Type, Length);
};

void TCigar::removeSoftClips(){
	for(auto it = _cigar.begin(); it!=_cigar.end(); ++it){
		if(it->type == 'S'){
			it = _cigar.erase(it);
		}
	}

	//update length
	_lengthSoftClippedLeft = 0;
	_lengthSoftClippedRight = 0;
};


}; //end namespace


