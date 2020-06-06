/*
 * TPSMCInput.h
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TPSMCINPUT_H_
#define GENOMETASKS_TPSMCINPUT_H_

#include "TGenome.h"
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
	GenotypeLikelihoods::TGenotypeData _prior;
	GenotypeLikelihoods::TGenotypePosteriorProbabilities _posterior;
	GenotypeLikelihoods::TBaseData _baseFreq;

	void _handleWindow();
public:
	TPSMCInput(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void createPSMCInput();
};

}; // end namespace

#endif /* GENOMETASKS_TPSMCINPUT_H_ */
