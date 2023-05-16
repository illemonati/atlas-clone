/*
 * TMajorMinor.h
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#ifndef TMAJORMINOR_H_
#define TMAJORMINOR_H_

#include <array>
#include <set>
#include <string>

#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "genometools/TGenotypeFrequencies.h"
#include "TGlfMultiReader.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TTask.h"
#include "coretools/Types/probability.h"

namespace PopulationTools {


//-----------------------------------------------
// TMajorMinorEstimator
//-----------------------------------------------


struct TMajorMinor {
	void run();
};

} // namespace PopulationTools

#endif /* TMAJORMINOR_H_ */
