/*
 * TestimateRecalibration.h
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TESTIMATERECALIBRATION_H_
#define GENOMETASKS_TESTIMATERECALIBRATION_H_

#include "TGenome.h"
#include "TRecalibrationEMEstimator.h"
#include "TTask.h"

namespace GenomeTasks{

//-----------------------------------------------------------
// TEstimateRecalibration_base
//-----------------------------------------------------------
class TEstimateRecalibration:public TGenome_windows{
protected:
	std::unique_ptr<GenotypeLikelihoods::RecalEstimator::TRecalibrationEMEstimator> recalObjectEM;
	BAM::TReadGroupMap* _readGroupMap;

	void _handleWindow();

public:
	TEstimateRecalibration(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	~TEstimateRecalibration();

	void estimateRecalibration();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_recal:public TTask{
public:
	TTask_recal(){
		_explanation = "Estimating error re-calibration parameters";
		_citations.insert("Kousathanas et al. (2017) Genetics");
	};

	void run(TParameters & Parameters, TLog* Logfile){
		TEstimateRecalibration estimator(Parameters, Logfile, _randomGenerator);
		estimator.estimateRecalibration();
	};
};


}; // end namespace

#endif /* GENOMETASKS_TESTIMATERECALIBRATION_H_ */
