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
#include "coretools/Types/TLogInt.h"
#include "coretools/Types/TPseudoInt.h"
#include "coretools/Types/probability.h"
#include "genometools/Genotypes/Base.h"
#include "genometools/Genotypes/BaseContext.h"

namespace BAM {

enum class End : size_t {min, from5=min, from3, max};
enum class Strand : size_t {min, Fwd=min, Rev, max};
enum class Mate : size_t {min, first=min, second, max};
enum class Flags: size_t {min, ReversedStrand = min, Paired, SecondMate, Aligned, SoftClipped, max};

inline std::string toString(Mate m) {
	constexpr coretools::TStrongArray<std::string_view, Mate> mates{{"Mate1", "Mate2"}};
	return std::string(mates[m]);
}

inline std::string toString(End e) {
	constexpr coretools::TStrongArray<std::string_view, End> ends{{"5", "3"}};
	return std::string(ends[e]);
}

inline std::string toString(Strand e) {
	constexpr coretools::TStrongArray<std::string_view, Strand> strands{{"Fwd", "Rev"}};
	return std::string(strands[e]);
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

	coretools::TLogInt fragmentLength = coretools::TLogInt::max();

	// distances are 0-based
	coretools::TPseudoInt distFrom5 = coretools::TPseudoInt::max();
	coretools::TPseudoInt distFrom3 = coretools::TPseudoInt::max();

	uint16_t readGroupID    = -1;
	uint16_t bamID          = -1;

	constexpr Mate mate() const noexcept {return static_cast<Mate>(get<Flags::SecondMate>());}
	constexpr End end() const noexcept {return distFrom3 < distFrom5 ? End::from3 : End::from5;}
	constexpr Strand strand() const noexcept {return get<Flags::ReversedStrand>() ? Strand::Rev : Strand::Fwd;}
	constexpr coretools::TPseudoInt dist(End E) const noexcept {return E==End::from5 ? distFrom5 : distFrom3;}
	constexpr genometools::BaseContext context() const {return genometools::baseContext(previousBase, base);}

	template<Flags F>
	constexpr bool get() const noexcept {return properties.get<F>();}
	template<Flags F>
	constexpr void set(bool B) noexcept {return properties.set<F>(B);}
};
}; // namespace BAM

#endif /* TBASE_H_ */
