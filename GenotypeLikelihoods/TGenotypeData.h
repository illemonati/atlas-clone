/*
 * TGenotypePrior.h
 *
 *  Created on: Nov 20, 2018
 *      Author: phaentu
 */

#ifndef TGENOTYPEDATA_H_
#define TGENOTYPEDATA_H_

#include <stddef.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <numeric>
#include <string>
#include <vector>

#include "TRandomGenerator.h"
#include "GenotypeTypes.h"
#include "probability.h"
#include "stringFunctions.h"
#include "TStrongArray.h"

namespace coretools { class TOutputFile; }

namespace GenotypeLikelihoods{

class Likelihood : public coretools::Probability{};

using TBaseProbability       = coretools::TStrongArray<coretools::Probability, genometools::Base, 4>;
using TBaseLikelikhoods      = coretools::TStrongArray<Likelihood, genometools::Base, 4>;
using TBaseData              = coretools::TStrongArray<double, genometools::Base, 4>;
using TBaseCounts            = coretools::TStrongArray<uint32_t, genometools::Base, 5>;
using TGenotypeLikelihoods   = coretools::TStrongArray<Likelihood, genometools::Genotype, 10>;
using TGenotypeProbabilities = coretools::TStrongArray<coretools::Probability, genometools::Genotype, 10>;
using TGenotypeData          = coretools::TStrongArray<double, genometools::Genotype, 10>;
}

#endif /* TGENOTYPEDATA_H_ */
