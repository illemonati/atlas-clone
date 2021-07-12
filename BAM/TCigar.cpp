/*
 * BamData.cpp
 *
 *  Created on: Jun 8, 2020
 *      Author: phaentu
 */

#include "TCigar.h"
#include "stringFunctions.h"

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
	_lengthSkipped = 0;
	_lengthSoftClippedLeft = 0;
	_lengthSoftClippedRight = 0;
	_addSoftClippedLeft = true;
	_softClippedAddedOnRight = false;
};

void TCigar::add(const char & Type, const uint32_t & Length){
	if(_softClippedAddedOnRight){
		throw "Cigar string contains entries past softclipping on right!";
	}
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
			_softClippedAddedOnRight = true;
		}
	} else if(Type == 'N'){
		_lengthSkipped += Length;
	} else if(Type != 'H' && Type != 'P'){
		//ignore
		std::string error = "Unknown CIGAR operation '";
		throw error + Type + "'!";
	}

	//add to vector
	_cigar.emplace_back(Type, Length);
};

void TCigar::removeSoftClips(){
	auto it = _cigar.begin();
	while(it!=_cigar.end()){
		if(it->type == 'S'){
			it = _cigar.erase(it);
		} else {
			++it;
		}
	}

	//update length
	_lengthSoftClippedLeft = 0;
	_lengthSoftClippedRight = 0;
};


}; //end namespace


