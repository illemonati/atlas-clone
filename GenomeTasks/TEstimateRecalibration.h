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

namespace BAM { class TReadGroupMap; }
namespace GenotypeLikelihoods { namespace SequencingError { class TRecalibrationEMEstimator; } }

namespace GenomeTasks{

//-----------------------------------------------------------
// TEstimateRecalibration_base
//-----------------------------------------------------------
class TEstimateRecalibration:public TGenome_windows{
protected:
	std::unique_ptr<GenotypeLikelihoods::SequencingError::TRecalibrationEMEstimator> recalObjectEM;
	BAM::TReadGroupMap* _readGroupMap;

	void _handleWindow() override;
	void _handleAlignment() override {}

public:
	TEstimateRecalibration();
	~TEstimateRecalibration();

	void estimateRecalibration();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_recal:public coretools::TTask{
public:
	TTask_recal(){
		_explanation = "Estimating error re-calibration parameters";
		_citations.insert("Kousathanas et al. (2017) Genetics");
	};

	void run(){
		using namespace coretools::instances;
		TEstimateRecalibration estimator;
		estimator.estimateRecalibration();
	};
};


}; // end namespace

#endif /* GENOMETASKS_TESTIMATERECALIBRATION_H_ */
