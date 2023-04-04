/*
 * TPMDEstimator.h
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TPMDESTIMATOR_H_
#define GENOMETASKS_TPMDESTIMATOR_H_

#include <memory>
#include <stdint.h>
#include <string>

#include "TGenome.h"
#include "PMD/TModels.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Main/TTask.h"

namespace BAM {
class TReadGroupMap;
}

namespace GenomeTasks {

//----------------------------------------
// TPMDEstimator.h
//----------------------------------------

class TPMDEstimator : public TGenome_parsed {
private:
	BAM::TReadGroupMap _readGroupMap;
	GenotypeLikelihoods::PMD::TModels* _pmd;
	void _handleAlignment();

public:
	TPMDEstimator();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TPMDESTIMATOR_H_ */
