/*
 * TGenotypePrior.h
 *
 *  Created on: Nov 20, 2018
 *      Author: phaentu
 */

#ifndef GENOTYPEDATA_H_
#define GENOTYPEDATA_H_

#include "coretools/Containers/TMassFunction.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Types/probability.h"

#include "genometools/GenotypeTypes.h"

namespace GenotypeLikelihoods{

using TBaseProbabilities     = coretools::TStrongMassFunction<coretools::Probability, genometools::Base, 4>;
using TBaseLikelihoods       = coretools::TStrongArray<coretools::Probability, genometools::Base, 4>;
using TBaseData              = coretools::TStrongArray<double, genometools::Base, 4>;
using TBaseCounts            = coretools::TStrongArray<uint32_t, genometools::Base, 5>;
using TGenotypeLikelihoods   = coretools::TStrongArray<coretools::Probability, genometools::Genotype, 10>;
using TGenotypeProbabilities = coretools::TStrongMassFunction<coretools::Probability, genometools::Genotype, 10>;
using TGenotypeData          = coretools::TStrongArray<double, genometools::Genotype, 10>;
using TBaseMassFunctions     = coretools::TStrongArray<TBaseProbabilities, genometools::Base, 4>;


} // namespace GenotypeLikelihoods

#endif /* TGENOTYPEDATA_H_ */
