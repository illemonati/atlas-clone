/*
 * TBase.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TSEQUENCEDBASE_H_
#define TSEQUENCEDBASE_H_

#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TBitSet.h"
#include <cstdint>
#include <iostream>

namespace BAM {

//---------------------------------------------------------------
// TSequencedBase
//---------------------------------------------------------------
// data container with minimal footprint, also used in recal
class TSequencedBase {
private:
	using bitset = coretools::TBitSet<3>;
	// flags: isReverseStrand, isSecondMate, isAligned
	bitset flags = 0; // initialized as 0,0,0
public:
	genometools::Base base           = genometools::N;
	genometools::BaseContext context = genometools::cNN;

	// original quality as in BAM file, but transformed to phredInt
	genometools::PhredIntProbability originalQuality_phredInt{0};
	// Quality after recalibration (used for filtering)
	genometools::PhredIntProbability recalibratedQualityAsPhredInt{0};

	// Do we need it if we also store fragment length?
	uint8_t mappingQuality  = 0;
	uint16_t fragmentLength = 0;
	uint16_t distFrom5Prime = -1; // zero based!
	uint16_t distFrom3Prime = -1; // zero based!
	uint16_t readGroupID    = -1;

	// set and get flags
	bool isReverseStrand() const noexcept { return flags.get<0>(); }
	bool isSecondMate() const noexcept { return flags.get<1>(); }
	bool isAligned() const noexcept { return flags.get<2>(); }

	void setReverseStrand(const bool status) noexcept  { flags.set<0>(status); }
	void setSecondMate(const bool status) noexcept  { flags.set<1>(status); }
	void setAligned(const bool status) noexcept  { flags.set<2>(status); }

	bool operator==(genometools::Base b) const noexcept { return base == b; }
	bool operator!=(genometools::Base b) const noexcept { return base != b; }
	void operator=(genometools::Base b) noexcept { base = b; }

	void print() const;
};

}; // namespace BAM

std::ostream &operator<<(std::ostream &os, BAM::TSequencedBase base);

#endif /* TBASE_H_ */
