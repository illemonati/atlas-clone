/*
 * TBase.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TSEQUENCEDBASE_H_
#define TSEQUENCEDBASE_H_

#include <cstdint>
#include <bitset>
#include <iostream>
#include "Types.h"

namespace BAM{

//---------------------------------------------------------------
//TSequencedBase
//---------------------------------------------------------------
//data container with minimal footprint, also used in recal
class TSequencedBase{
private:
	//flags: isReverseStrand, isSecondMate, isAligned
	std::bitset<3> flags; //initialized as 0,0,0

public:
	Base base;
	PhredIntErrorRate originalQuality_phredInt; //original quality as in BAM file, but transformed to phredInt
	PhredIntErrorRate recalibratedQualityAsPhredInt; //Quality after recalibration (used for filtering)
	uint16_t distFrom5Prime; //zero based!
	uint16_t distFrom3Prime; //zero based!	Do we need it if we also store fragment length?
	uint16_t readGroupID;
	uint16_t fragmentLength;
	uint8_t mappingQuality;
	BaseContext context;

	TSequencedBase();

	//set and get flags
	bool isReverseStrand() const { return flags[0]; };
	bool isSecondMate() const { return flags[1]; };
	bool isAligned() const { return flags[2]; };

	void setReverseStrand(const bool status){ flags[0] = status; };
	void setSecondMate(const bool status){ flags[1] = status; };
	void setAligned(const bool status){ flags[2] = status; };

	bool operator==(const Base & b) const { return base == b; };
	bool operator!=(const Base & b) const { return base != b; };

	void operator=(const Base & b){ base = b; };

	void print() const;
};

}; //end namespace

std::ostream& operator<<(std::ostream& os, const BAM::TSequencedBase & base);

#endif /* TBASE_H_ */
