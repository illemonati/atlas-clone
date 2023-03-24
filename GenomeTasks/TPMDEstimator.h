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
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TLog.h"
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
	size_t _tableSize;
	BAM::TReadGroupMap _readGroupMap;
	GenotypeLikelihoods::TPostMortemDamage::PMDTables _pmdTables;
	void _handleAlignment();
	void _write(bool normalized);

public:
	TPMDEstimator();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TPMDESTIMATOR_H_ */
