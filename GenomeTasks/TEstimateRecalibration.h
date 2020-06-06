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

namespace GenomeTasks{

//-----------------------------------------------------------
// TEstimateRecalibration_base
//-----------------------------------------------------------
class TEstimateRecalibration_base:public TGenome_windows{
private:
	std::unique_ptr<GenotypeLikelihoods::TRecalibrationEMEstimator> recalObjectEM;
	std::unique_ptr<BAM::TReadGroupMap> _readGroupMap;
	GenotypeLikelihoods::TBaseData _baseFreq;

	void _handleWindow();

public:
	TEstimateRecalibration_base(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
};


//-----------------------------------------------------------
// TEstimateRecalibration
//-----------------------------------------------------------
class TEstimateRecalibration:public TEstimateRecalibration_base{
private:
	bool _writeTmpTables;

public:
	TEstimateRecalibration(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void estimateRecalibration();
};

//-----------------------------------------------------------
// TEstimateRecalibrationLL
//-----------------------------------------------------------
class TEstimateRecalibrationLL:public TEstimateRecalibration_base{
public:
	TEstimateRecalibrationLL(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void calculateRecalibrationLL();
};

}; // end namespace

#endif /* GENOMETASKS_TESTIMATERECALIBRATION_H_ */
