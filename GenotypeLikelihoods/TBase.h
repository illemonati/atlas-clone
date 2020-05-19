/*
 * TBase.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TBASE_H_
#define TBASE_H_

#include "TGenotypeMap.h"
#include <bitset>


//---------------------------------------------------------------
//TBase
//---------------------------------------------------------------
//data container with minimal footprint, also used in recal
class TBase{
private:
	//flags: isReverseStrand, isSecondMate, isAligned
	std::bitset<3> flags; //initialized as 0,0,0

public:
	Base base;
	uint8_t originalQuality_phredInt; //original quality as in BAM file, but transformed to phredInt
	uint8_t recalibratedQualityAsPhredInt; //Quality after recalibration (used for filtering)
	uint16_t distFrom5Prime; //zero based!
	uint16_t distFrom3Prime; //zero based!	Do we need it if we also store fragment length?
	uint16_t readGroup;
	uint16_t fragmentLength;
	uint8_t mappingQuality;
	BaseContext context;

	TBase(){
		base = N;
		originalQuality_phredInt = 0;
		recalibratedQualityAsPhredInt = 0;
		distFrom5Prime = -1;
		distFrom3Prime = -1;
		readGroup = -1;
		context = cNN;
		mappingQuality = 0;
		fragmentLength = 0;
	};

	//set and get flags
	bool isReverseStrand() const { return flags[0]; };
	bool isSecondMate() const { return flags[1]; };
	bool isAligned() const { return flags[2]; };

	void setReverseStrand(const bool status){ flags[0] = status; };
	void setSecondMate(const bool status){ flags[1] = status; };
	void setAligned(const bool status){ flags[2] = status; };
};

#endif /* TBASE_H_ */
