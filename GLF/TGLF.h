/*
 * TGLF.h
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#ifndef TGLF_H_
#define TGLF_H_

#include "coretools/Containers/TDualStrongArray.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "genometools/VCF/TVcfWriter.h"

namespace GLF {

using TGLFLikelihoods = coretools::TDualStrongArray<genometools::HighPrecisionPhredIntProbability, genometools::Base,
                                                    genometools::Genotype, 4, 10, genometools::Ploidy>;

constexpr std::string_view version() noexcept { return "GLF2"; }

//--------------------------------------------
// Tasks
//--------------------------------------------
struct TGLFPrinter {
	void run();
};

struct TGLFIndexer {
	void run();
};

}; // namespace GLF

#endif /* TGLF_H_ */
