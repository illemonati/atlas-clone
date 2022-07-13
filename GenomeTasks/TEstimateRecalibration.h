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
#include "TTask.h"
#include "TRecalibrationEMEstimator.h"

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
	std::unique_ptr<GenotypeLikelihoods::SequencingError::TRecalibrationEMEstimator> recalObjectEM;
	std::unique_ptr<BAM::TReadGroupMap> _readGroupMap;
	bool _onlyLL = false;

	void _handleWindow() override;
	void _handleAlignment() override {}

public:
	TEstimateRecalibration();

	void estimateRecalibration();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_recal : public coretools::TTask {
public:
	TTask_recal() {
		_explanation = "Estimating error re-calibration parameters";
		_citations.insert("Kousathanas et al. (2017) Genetics");
	};

	void run() {
		using namespace coretools::instances;
		TEstimateRecalibration estimator;
		estimator.estimateRecalibration();
	};
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TESTIMATERECALIBRATION_H_ */
