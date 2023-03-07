/*
 * TestimateRecalibration.h
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TESTIMATERECALIBRATION_H_
#define GENOMETASKS_TESTIMATERECALIBRATION_H_

#include <memory>
#include <set>
#include <string>

#include "TGenome.h"
#include "coretools/Main/TTask.h"
#include "SequencingError/TRecalEstimator.h"

namespace BAM {
class TReadGroupMap;
}
namespace GenotypeLikelihoods {
namespace SequencingError {
class TRecalibrationEMEstimator;
}
} // namespace GenotypeLikelihoods

namespace GenomeTasks {

//-----------------------------------------------------------
// TEstimateRecalibration_base
//-----------------------------------------------------------
class TEstimateRecalibration final : public TGenome_windows {
private:
	BAM::TReadGroupMap _readGroupMap;
	GenotypeLikelihoods::SequencingError::TRecalibrationEMEstimator recal;
	bool _onlyLL = false;

	void _handleWindow() override;
	void _handleAlignment() override {}

public:
	TEstimateRecalibration();

	void run();
};

} // namespace GenomeTasks

#endif /* GENOMETASKS_TESTIMATERECALIBRATION_H_ */
