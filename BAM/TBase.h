/*
 * TBase.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TBASE_H_
#define TBASE_H_

#include <bitset>
#include <iostream>

enum Base : uint8_t {A=0, C, G, T, N};

namespace BAM{

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
	uint16_t readGroupID;
	uint16_t fragmentLength;
	uint8_t mappingQuality;
	Base context;

	TBase();

	//set and get flags
	bool isReverseStrand() const { return flags[0]; };
	bool isSecondMate() const { return flags[1]; };
	bool isAligned() const { return flags[2]; };

	void setReverseStrand(const bool status){ flags[0] = status; };
	void setSecondMate(const bool status){ flags[1] = status; };
	void setAligned(const bool status){ flags[2] = status; };

	bool operator==(const Base & b){ return base == b; };
	bool operator!=(const Base & b){ return base != b; };

	void operator=(const Base & b){ base = b; };

	void print();
};

}; //end namespace


std::ostream& operator<<(std::ostream& os, const BAM::TBase & base);
std::ostream& operator<<(std::ostream& os, const Base & base);


#endif /* TBASE_H_ */
