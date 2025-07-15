/*
 * TMajorMinor.h
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#ifndef TMAJORMINOR_H_
#define TMAJORMINOR_H_

#include "coretools/Main/TError.h"
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
std::array<genometools::AllelicCombination, 3> useAllelicCombinationsThatContain(genometools::Base base);;
coretools::Log10Probability LLFixedAllele(coretools::TConstView<genometools::TGLFEntry> data, genometools::Base major);

}

struct TMajorMinor {
	void run();
};

} // namespace PopulationTools

#endif /* TMAJORMINOR_H_ */
