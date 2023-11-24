/*
 * TBase.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TSEQUENCEDBASE_H_
#define TSEQUENCEDBASE_H_

#include "coretools/Containers/TBitSet.h"
#include "coretools/Containers/TStrongArray.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"

namespace BAM {

enum class Mate : size_t {min, first=min, second, max};
inline std::string toString(Mate m) {
	constexpr coretools::TStrongArray<std::string_view, Mate> mates{{"Mate1", "Mate2"}};
	return std::string(mates[m]);
}

//---------------------------------------------------------------
// TSequencedBase
//---------------------------------------------------------------
// data container with minimal footprint, also used in recal
class TSequencedBase {
public:
	genometools::Base base         = genometools::Base::N;
	genometools::Base previousBase = genometools::Base::N;

	// original quality as in BAM file, but transformed to phredInt
	genometools::PhredIntProbability originalQuality_phredInt{0};
	// Quality after recalibration (used for filtering)
	genometools::PhredIntProbability recalibratedQualityAsPhredInt{0};

	coretools::TBitSet<3> flags{0}; // initialized as 0,0,0

	genometools::PhredIntProbability mappingQuality{0};

	// Do we need it if we also store fragment length?
	uint16_t fragmentLength = 0;
	uint16_t distFrom5Prime = -1; // zero based!
	uint16_t distFrom3Prime = -1; // zero based!
	uint16_t readGroupID    = -1;

	// set and get flags
	bool isReverseStrand() const noexcept { return flags.get<0>(); }
	bool isSecondMate() const noexcept { return flags.get<1>(); }
	Mate mate() const noexcept {return static_cast<Mate>(flags.get<1>());}
	bool isAligned() const noexcept { return flags.get<2>(); }

	void setReverseStrand(const bool status) noexcept  { flags.set<0>(status); }
	void setSecondMate(const bool status) noexcept  { flags.set<1>(status); }
	void setAligned(const bool status) noexcept  { flags.set<2>(status); }

	bool operator==(genometools::Base b) const noexcept { return base == b; }
	bool operator!=(genometools::Base b) const noexcept { return base != b; }
	void operator=(genometools::Base b) noexcept { base = b; }
	operator bool() const noexcept {return base != genometools::Base::N;}

	constexpr genometools::BaseContext context() const {return genometools::baseContext(previousBase, base);}

	void print() const;
};

}; // namespace BAM

std::ostream &operator<<(std::ostream &os, BAM::TSequencedBase base);

#endif /* TBASE_H_ */
