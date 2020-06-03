/*
 * TPMDSCalculator.h
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TPMDSCALCULATOR_H_
#define GENOMETASKS_TPMDSCALCULATOR_H_

#include "TGenome.h"

//----------------------------------------------
// TPMDSCalculator
//----------------------------------------------
class TPMDSCalculator:public TGenome_parsed{
private:
	double _pi, _minPMDS, _maxPMDS;
	BAM::TOutputBamFile _outBam;

	void _handleAlignment();

public:
	TPMDSCalculator(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void estimatePMDS();
};



#endif /* GENOMETASKS_TPMDSCALCULATOR_H_ */
