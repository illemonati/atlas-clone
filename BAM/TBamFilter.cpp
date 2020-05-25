/*
 * TBamFilter.cpp
 *
 *  Created on: May 24, 2020
 *      Author: phaentu
 */


#include "TBamFilter.h"

namespace BAM{

//-----------------------------------------------------
//TBamFileFilter_base
//-----------------------------------------------------
TBamFileFilter_base::TBamFileFilter_base(){
	_counter = 0;
	_keep = true;
	_updateBlacklist = false;
	_blacklist = nullptr;
};

void TBamFileFilter_base::_filterOut(const std::string & alignmentName, const bool & isReverseStrand){
	++_counter;
	if(_updateBlacklist){
		_blacklist->addToBlacklist(alignmentName, isReverseStrand, _errorString);
	}
};

void TBamFileFilter_base::keep(){
	_keep = true;
};

void TBamFileFilter_base::setBlacklist(const TAlignmentBlacklist* Blacklist){
	_blacklist = Blacklist;
	_updateBlacklist = true;
};

//-----------------------------------------------------
//TBamFileFilterBool
//-----------------------------------------------------
void TBamFileFilterBool::filter(const std::string ErrorString){
	_keep = false;
	_errorString = ErrorString;
};

bool TBamFileFilterBool::pass(const bool state, const std::string & alignmentName, const bool & isReverseStrand){
	if(state && !_keep){
		_filterOut(alignmentName, isReverseStrand);
		return false;
	}
	return true;
};

//-----------------------------------------------------
//TBamFileFilterRange
//-----------------------------------------------------
TBamFileFilterRange::TBamFileFilterRange(){
	_min = 0;
	_max = UINT32_MAX;
};

void TBamFileFilterRange::filter(const uint32_t Min, const uint32_t Max, const std::string ErrorString){
	_keep = false;
	_min = Min;
	_max = Max;
	_errorString = ErrorString;
};

bool TBamFileFilterRange::pass(const uint32_t value, const std::string & alignmentName, const bool & isReverseStrand){
	if(!_keep && (value < _min || value > _max)){
		_filterOut(alignmentName, isReverseStrand);
	}
	return true;
};


}; //end namespace
