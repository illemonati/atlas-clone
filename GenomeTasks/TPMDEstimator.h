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
#include "coretools/Main/TLog.h"
#include "TPMDTables.h"
#include "coretools/Main/TParameters.h"
#include "TPostMortemDamage.h"
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
	uint16_t _maxLengthForInference;
	std::unique_ptr<BAM::TReadGroupMap> _readGroupMap;
	GenotypeLikelihoods::TPMDTables _pmdTables;
	GenotypeLikelihoods::TPMDEstimationParameters _estimationParameters;

	void _handleAlignment();

public:
	TPMDEstimator();
	void estimatePMD();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_estimatePMD : public coretools::TTask {
public:
	TTask_estimatePMD() { _explanation = "Estimating Post-Mortem Damage (PMD) patterns"; };

	void run() {
		using namespace coretools::instances;
		TPMDEstimator estimator;
		estimator.estimatePMD();
	};
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TPMDESTIMATOR_H_ */
