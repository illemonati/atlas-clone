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
class TEstimateRecalibration_base:public TGenome_windows{
protected:
	std::unique_ptr<GenotypeLikelihoods::TRecalibrationEMEstimator> recalObjectEM;
	BAM::TReadGroupMap* _readGroupMap;
	GenotypeLikelihoods::TBaseData _baseFreq;

	void _handleWindow();

public:
	TEstimateRecalibration_base(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	virtual ~TEstimateRecalibration_base();
};


//-----------------------------------------------------------
// TEstimateRecalibration
//-----------------------------------------------------------
class TEstimateRecalibration:public TEstimateRecalibration_base{
private:
	bool _writeTmpTables;

public:
	TEstimateRecalibration(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void estimateRecalibration();
};

//-----------------------------------------------------------
// TEstimateRecalibrationLL
//-----------------------------------------------------------
class TEstimateRecalibrationLL:public TEstimateRecalibration_base{
public:
	TEstimateRecalibrationLL(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void calculateRecalibrationLL();
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

class TTask_recalLL:public TTask{
public:
	TTask_recalLL(){
		_explanation = "Calculating LL for error re-calibration function";
		_citations.insert("Kousathanas et al. (2017) Genetics");
	};

	void run(TParameters & Parameters, TLog* Logfile){
		TEstimateRecalibrationLL estimator(Parameters, Logfile, _randomGenerator);
		estimator.calculateRecalibrationLL();
	};
};



}; // end namespace

#endif /* GENOMETASKS_TESTIMATERECALIBRATION_H_ */
