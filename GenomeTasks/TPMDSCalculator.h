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
		TPMDSCalculator calculator;
		calculator.calculatePMDS();
	};
};

}; // end namespace


#endif /* GENOMETASKS_TPMDSCALCULATOR_H_ */
