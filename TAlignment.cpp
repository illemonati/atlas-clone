/*
 * TAlignment.cpp
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */

#include "TAlignment.h"


void TAlignment::clear(){
	empty = true;
	parsed = false;
	recalibrated = false;
	changed = false;
	passedFilters = false;
	quality = qualityOriginal;
}

void TAlignment::fill(BamTools::BamAlignment & bamAlignment, int ReadGroupId){
	//clear
	clear();

	//add basic info
	chrNumber = bamAlignment.RefID;
	position = bamAlignment.Position;
	isReverseStrand = bamAlignment.IsReverseStrand();
	isProperPair = bamAlignment.IsProperPair();
	readGroupId = ReadGroupId;
}

void TAlignment::setFiltersPassed(bool passed){
	passedFilters = passed;
};


