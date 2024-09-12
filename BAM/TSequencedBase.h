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
#include "coretools/Types/probability.h"
#include "genometools/Genotypes/Base.h"
#include "genometools/Genotypes/BaseContext.h"

namespace BAM {

enum class End : size_t {min, from5=min, from3, max};
enum class Mate : size_t {min, first=min, second, max};
enum class Flags: size_t {min, ReversedStrand = min, SecondMate, Aligned, SoftClipped, max};

inline std::string toString(Mate m) {
	constexpr coretools::TStrongArray<std::string_view, Mate> mates{{"Mate1", "Mate2"}};
	return std::string(mates[m]);
}

//---------------------------------------------------------------
// TSequencedBase
//---------------------------------------------------------------
// data container with minimal footprint, also used in recal
struct TSequencedBase {
	genometools::Base base         = genometools::Base::N;
	genometools::Base previousBase = genometools::Base::N;

	// original quality as in BAM file, but transformed to phredInt
	coretools::PhredInt originalQuality = coretools::PhredInt::highest();
	// Quality after recalibration (used for filtering)
	coretools::PhredInt recalQuality = coretools::PhredInt::highest();

	coretools::PhredInt mappingQuality = coretools::PhredInt::highest();

	coretools::TStrongBitSet<Flags> properties{0}; // initialized as 0,0,0,0

	uint16_t fragmentLength = 0;
	uint16_t distFrom5      = 0; // zero based!
	uint16_t distFrom3      = 0; // zero based!

	uint16_t readGroupID    = -1;
	uint16_t bamID          = -1;

	constexpr Mate mate() const noexcept {return static_cast<Mate>(get<Flags::SecondMate>());}
	constexpr End end() const noexcept {return distFrom5 < distFrom3 ? End::from5 : End::from3;}
	constexpr uint16_t dist(End E) const noexcept {return E==End::from5 ? distFrom5 : distFrom3;}
	constexpr genometools::BaseContext context() const {return genometools::baseContext(previousBase, base);}

	template<Flags F>
	constexpr bool get() const noexcept {return properties.get<F>();}
	template<Flags F>
	constexpr void set(bool B) noexcept {return properties.set<F>(B);}
};
}; // namespace BAM

#endif /* TBASE_H_ */
