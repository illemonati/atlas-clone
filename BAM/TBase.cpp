/*
 * TBase.cpp
 *
 *  Created on: Jul 31, 2020
 *      Author: wegmannd
 */

#include "TBase.h"

namespace BAM{

TBase::TBase(){
	base = N;
	originalQuality_phredInt = 0;
	recalibratedQualityAsPhredInt = 0;
	distFrom5Prime = -1;
	distFrom3Prime = -1;
	readGroupID = -1;
	context = N;
	mappingQuality = 0;
	fragmentLength = 0;
};

std::ostream& operator<<(std::ostream& os, const TBase & base){
	os << base.base;
	return os;
};

std::ostream& operator<<(std::ostream& os, const Base & base){
	if(base == A){
		os << "A";
	} else if(base == C){
		os << 'C';
	} else if(base == G){
		os << 'G';
	} else if(base == T){
		os << 'T';
	} else {
		os << 'N';
	}
	return os;
};

}; //end namespace

