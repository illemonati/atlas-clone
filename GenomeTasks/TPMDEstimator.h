/*
 * TPMDEstimator.h
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TPMDESTIMATOR_H_
#define GENOMETASKS_TPMDESTIMATOR_H_

#include "TBamTraverser.h"
#include "oldPMD/TModels.h"

namespace BAM {
class TReadGroupMap;
}

namespace GenomeTasks {

//----------------------------------------
// TPMDEstimator.h
//----------------------------------------

class TPMDEstimator final : public TBamTraverser<true> {
private:
	BAM::TReadGroupMap _readGroupMap;
	GenotypeLikelihoods::oldPMD::TModels _pmd;
	void _handleAlignment(BAM::TAlignment& alignment) override;

public:
	TPMDEstimator();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TPMDESTIMATOR_H_ */
