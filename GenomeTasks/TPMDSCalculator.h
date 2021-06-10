/*
 * TPMDSCalculator.h
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TPMDSCALCULATOR_H_
#define GENOMETASKS_TPMDSCALCULATOR_H_

#include "TGenome.h"
#include "TTask.h"

namespace GenomeTasks{

//----------------------------------------------
// TPMDSCalculator
//----------------------------------------------
class TPMDSCalculator:public TGenome_parsed{
private:
	double _pi, _minPMDS, _maxPMDS;
	BAM::TOutputBamFile _outBam;

	double _calculatePMDS();
	void _handleAlignment();

public:
	TPMDSCalculator(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void calculatePMDS();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_PMDS:public coretools::TTask{
public:
	TTask_PMDS(){
		_explanation = "Filtering for ancient reads using PMDS";
		_citations.insert("Skoglund et al. (2014) PNAS");
	};

	void run(TParameters & Parameters, TLog* Logfile){
		TPMDSCalculator calculator(Parameters, Logfile, _randomGenerator);
		calculator.calculatePMDS();
	};
};

}; // end namespace


#endif /* GENOMETASKS_TPMDSCALCULATOR_H_ */
