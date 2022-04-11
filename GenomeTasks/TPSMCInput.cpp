/*
 * TPSMCInput.cpp
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#include "TPSMCInput.h"

#include <math.h>
#include <exception>
#include <ostream>

#include "TGenotypeData.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TSite.h"
#include "TThetaEstimator.h"
#include "TWindow.h"
#include "probability.h"

namespace GenomeTasks{

//----------------------------------------
// TPSMCInput
//----------------------------------------
TPSMCInput::TPSMCInput(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator):TGenome_windows(Parameters, Logfile, RandomGenerator){
	_theta = Parameters.getParameterWithDefault<double>("theta", 0.001);
	_logfile->list("Using theta = ", _theta, ". (parameter 'theta')");

	_thetaEstimator = std::make_unique<GenotypeLikelihoods::TThetaEstimator>(_logfile, _randomGenerator);
	_thetaEstimator->setTheta(_theta);

	_confidence = Parameters.getParameterWithDefault<double>("confidence", 0.99);
	_logfile->list("Calling heterozygosity state with confidence > ", _confidence, ". (parameter 'confidence')");
	_logConfidence = log(_confidence);
	_logConfidenceHet = log(1.0 - _confidence);


	_blockSize = Parameters.getParameterWithDefault<int>("block", 100);
	//make sure window size is a multiple of block length!
	if(_windowSize % _blockSize != 0) throw "Window size is not a multiple of block size!";
	_nBlocks = _window.size() / _blockSize;

	//open output file
	std::string outputFileName = _outputName + ".psmcfa";
	_logfile->list("Writing PSMC input file to '" + outputFileName + "'.");
	_out.open(outputFileName.c_str());
	if(!_out) throw "Failed to open output file '" + outputFileName + "'!";
	_nCharOnLine = 0;
};

void TPSMCInput::_handleWindow(){
	_logfile->listFlushTime("Estimating heterozygosity status ...");

	//calc prior probabilities on Genotypes
	_thetaEstimator->fillPGenotype(_prior);

	//estimating base frequencies
	_thetaEstimator->setBaseFreq( _window.estimateBaseFrequencies() );

	//call heterozygosity in blocks
	for(uint32_t b=0; b<_nBlocks; ++b){
		uint32_t blockStart = b * _blockSize;
		coretools::LogProbability logPHomo;

		for(uint32_t i=0; i<_blockSize; ++i){
			if(!_window[blockStart + i].empty()){
				_genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(_window[blockStart + 1], _genoLik);
				_posterior = GenotypeLikelihoods::posterior(_genoLik, _prior);
				logPHomo += coretools::LogProbability(GenotypeLikelihoods::homozygous(_posterior));
			}
		}

		//check if we are heterozygous
		if(logPHomo > _logConfidence){
			_out << 'T';
		} else if(logPHomo < _logConfidenceHet){
			_out << 'K';
		} else {
			_out << 'N';
		}

		//do we add a new line?
		if(_nCharOnLine == 59){
			_nCharOnLine = 0;
			_out << '\n';
		} else ++_nCharOnLine;
	}

	_logfile->doneTime();
};

void TPSMCInput::createPSMCInput(){
	_traverseBAMWindows();

	//close output file
	if(_nCharOnLine > 0) _out << '\n';
	_out.close();
};

}; // end namespace

