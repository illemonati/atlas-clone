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
	coretools::Probability _pi;
	double _minPMDS, _maxPMDS;
	BAM::TOutputBamFile _outBam;

	double _calculatePMDS();
	void _handleAlignment();

public:
	TPMDSCalculator(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
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

	void run(){
		using namespace coretools::instances;
		TPMDSCalculator calculator(parameters(), &logfile(), &randomGenerator());
		calculator.calculatePMDS();
	};
};

}; // end namespace


#endif /* GENOMETASKS_TPMDSCALCULATOR_H_ */
