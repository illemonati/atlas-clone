/*
 * TPMDEstimator.h
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TPMDESTIMATOR_H_
#define GENOMETASKS_TPMDESTIMATOR_H_

#include "TGenome.h"
#include "TPostMortemDamage.h"

//----------------------------------------
// TPMDEstimator.h
//----------------------------------------

class TPMDEstimator:public TGenome_parsed{
private:
	uint16_t _maxLengthForInference;
	BAM::TReadGroupMap* _readGroupMap;
	GenotypeLikelihoods::TPMDTables _pmdTables;
	bool _estimateExponential;
	double _eps;
	uint16_t _numNRIterations;

	void _handleAlignment();

public:
	TPMDEstimator(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	~TPMDEstimator();
	void estimatePMD();
};


#endif /* GENOMETASKS_TPMDESTIMATOR_H_ */
