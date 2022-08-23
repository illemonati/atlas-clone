/*
 * TBase.cpp
 *
 *  Created on: Jul 31, 2020
 *      Author: wegmannd
 */

#include "TSequencedBase.h"

namespace BAM {

void TSequencedBase::print() const {
	std::cout << "base: " << base << '\n'
		  << "originalQuality_phredInt: " << originalQuality_phredInt << '\n'
		  << "recalibratedQualityAsPhredInt: " << recalibratedQualityAsPhredInt << '\n'
		  << "distFrom5Prime: " << distFrom5Prime << '\n'
		  << "distFrom3Prime: " << distFrom3Prime << '\n'
		  << "readGroupID: " << readGroupID << '\n'
			  << "context: " << context() << '\n'
		  << "mappingQuality: " << mappingQuality << '\n'
		  << "fragmentLength: " << fragmentLength << std::endl;
}

}; // namespace BAM

std::ostream &operator<<(std::ostream &os, BAM::TSequencedBase base) {
	os << base.base;
	return os;
}
