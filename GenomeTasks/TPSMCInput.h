/*
 * TPSMCInput.h
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TPSMCINPUT_H_
#define GENOMETASKS_TPSMCINPUT_H_

#include <stdint.h>
#include <iosfwd>
#include <memory>
#include <string>

#include "TGenome.h"
#include "TGenotypeData.h"
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TTask.h"
#include "TThetaEstimator.h"

namespace GenomeTasks{

//----------------------------------------
// TPSMCInput
//----------------------------------------
class TPSMCInput:public TGenome_windows{
private:
	double _theta;
	double _confidence, _logConfidence, _logConfidenceHet;

	uint32_t _blockSize;
	uint32_t _nBlocks;
	std::ofstream _out;
	uint32_t _nCharOnLine;
	std::unique_ptr<GenotypeLikelihoods::TThetaEstimator> _thetaEstimator;

	//tmp variables
	GenotypeLikelihoods::TGenotypeProbabilities _prior;
	GenotypeLikelihoods::TGenotypeProbabilities _posterior;

	void _handleWindow();
public:
	TPSMCInput(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void createPSMCInput();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_PSMC:public coretools::TTask{
public:
	TTask_PSMC(){ _explanation = "Generating a PSMC Input file probabilistically"; };

	void run(){
		using namespace coretools::instances;
		TPSMCInput psmc(parameters(), &logfile(), &randomGenerator());
		psmc.createPSMCInput();
	};
};


}; // end namespace

#endif /* GENOMETASKS_TPSMCINPUT_H_ */
