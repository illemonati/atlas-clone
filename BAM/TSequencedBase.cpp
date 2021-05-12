/*
 * TBase.cpp
 *
 *  Created on: Jul 31, 2020
 *      Author: wegmannd
 */

#include "TSequencedBase.h"

namespace BAM{

TSequencedBase::TSequencedBase(){
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

void TSequencedBase::print() const{
	std::cout << "base: " << base << std::endl;
	std::cout << "originalQuality_phredInt: " << originalQuality_phredInt << std::endl;
	std::cout << "recalibratedQualityAsPhredInt: " << recalibratedQualityAsPhredInt << std::endl;
	std::cout << "distFrom5Prime: " << distFrom5Prime << std::endl;
	std::cout << "distFrom3Prime: " << distFrom3Prime << std::endl;
	std::cout << "readGroupID: " << readGroupID << std::endl;
	std::cout << "context: " << context << std::endl;
	std::cout << "mappingQuality: " << mappingQuality << std::endl;
	std::cout << "fragmentLength: " << fragmentLength << std::endl;
};

}; //end namespace

std::ostream& operator<<(std::ostream& os, const BAM::TSequencedBase & base){
	os << base.base;
	return os;
};
