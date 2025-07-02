/*
 * TMajorMinor.h
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#ifndef TMAJORMINOR_H_
#define TMAJORMINOR_H_

#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"
#include "coretools/enum.h"

#include "genometools/GLF/GLF.h"
#include "genometools/Genotypes/AllelicCombination.h"
#include "genometools/Genotypes/Base.h"

namespace PopulationTools {

namespace MajorMinor {
struct TMMData {
	bool pass{false};
	coretools::Probability MAF;
	coretools::PhredInt variantQuality;
	genometools::Base major;
	genometools::Base minor;
};

constexpr TMMData failedTMMData = {false, coretools::Probability::lowest(), coretools::PhredInt::lowest(), genometools::Base::N, genometools::Base::N};

constexpr auto useAllelicCombinationsThatContain(genometools::Base base) {
	using genometools::Base;
	assert(base != Base::N);
	using AC = genometools::AllelicCombination;
	switch (base) {
	case Base::A: return std::array{AC::AC, AC::AG, AC::AT};
	case Base::C: return std::array{AC::AC, AC::CG, AC::CT};
	case Base::G: return std::array{AC::AG, AC::CG, AC::GT};
	default: return std::array{AC::AT, AC::CT, AC::GT};
	}
};

template<typename Container> genometools::AllelicCombination chooseBestAllelicCombination(const Container &acd) {
	return coretools::instances::randomGenerator().sampleIndexOfMaxima<Container, genometools::AllelicCombination>(acd);
};

coretools::Log10Probability LLFixedAllele(coretools::TConstView<genometools::TGLFEntry> data, genometools::Base major);

}

struct TMajorMinor {
	void run();
};

} // namespace PopulationTools

#endif /* TMAJORMINOR_H_ */
