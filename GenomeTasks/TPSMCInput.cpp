/*
 * TPSMCInput.cpp
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#include "TPSMCInput.h"

namespace GenomeTasks{

//----------------------------------------
// TPSMCInput
//----------------------------------------
TPSMCInput::TPSMCInput(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_windows(Params, Logfile, RandomGenerator){
	_theta = Params.getParameterDoubleWithDefault("theta", 0.001);
	_logfile->list("Using theta = " + toString(_theta) + ". (parameter 'theta')");

	_thetaEstimator = std::make_unique<TThetaEstimator>(_logfile, _randomGenerator);
	_thetaEstimator->setTheta(_theta);

	_confidence = Params.getParameterDoubleWithDefault("confidence", 0.99);
	_logfile->list("Calling heterozygosity state with confidence > " + toString(_confidence) + ". (parameter 'confidence')");
	_logConfidence = log(_confidence);
	_logConfidenceHet = log(1.0 - _confidence);


	_blockSize = Params.getParameterIntWithDefault("block", 100);
	//make sure window size is a multiple of block length!
	if(_windowSize % _blockSize != 0) throw "Window size is not a multiple of block size!";
	_nBlocks = _window.length / _blockSize;

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
	_window.estimateBaseFrequencies(_baseFreq);
	_thetaEstimator->setBaseFreq(_baseFreq);

	//call heterozygosity in blocks
	for(int b=0; b<_nBlocks; ++b){
		uint32_t blockStart = b * _blockSize;
		double logPHomo = 0.0;

		for(int i=0; i<_blockSize; ++i){
			if(_window[blockStart + i].hasData){
				_genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(_window[blockStart + 1].bases, _genoLik);
				_posterior.fill(_genoLik, _prior);
				logPHomo += log(_posterior.probHomozygous());
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

