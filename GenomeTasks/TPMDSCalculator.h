/*
 * TPMDSCalculator.h
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TPMDSCALCULATOR_H_
#define GENOMETASKS_TPMDSCALCULATOR_H_

#include <set>
#include <string>

#include "TBamFile.h"
#include "TGenome.h"
#include "coretools/Main/TTask.h"
#include "coretools/Types/probability.h"

namespace GenomeTasks{

//----------------------------------------------
// TPMDSCalculator
//----------------------------------------------
class TPMDSCalculator:public TGenome_parsed{
private:
	coretools::Probability _pi;
	double _minPMDS, _maxPMDS;
	BAM::TOutputBamFile _outBam;

	double _calculatePMDS();
	void _handleAlignment() override;

public:
	TPMDSCalculator();
	void run();
};

}; // end namespace


#endif /* GENOMETASKS_TPMDSCALCULATOR_H_ */
