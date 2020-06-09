/*
 * TPMDEstimator.cpp
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#include "TPMDEstimator.h"

namespace GenomeTasks{

//----------------------------------------
// TPMDEstimator.h
//----------------------------------------
TPMDEstimator::TPMDEstimator(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_parsed(Parameters, Logfile, RandomGenerator){
	//make sure it has a reference
	_openReference(true);

	//prepare maps
	_readGroupMap = new BAM::TReadGroupMap(&_bamFile._readGroups, Parameters.getParameterString("poolReadGroups", false), _logfile);

	//prepare PMD table
	_maxLengthForInference = Parameters.getParameterIntWithDefault("length", 50);
	_logfile->list("Estimating PMD at the first " + toString(_maxLengthForInference) + " positions.");

	if(Parameters.parameterExists("onlyEmpiric")){
		_logfile->list("Not fitting exponential models. (parameter 'onlyEmpiric')");
		_estimateExponential = false;
	} else {
		_logfile->startIndent("Will estimate exponential PMD models: (use 'onlyEmpiric' to turn off)");
		_eps = Parameters.getParameterDoubleWithDefault("eps", 0.001);
		_logfile->list("Will conider the Newton-Raphson algorithm to have converged if the likelihood difference < " + toString(_eps) + ". (parameter 'eps')");
		_numNRIterations = Parameters.getParameterIntWithDefault("numNRIterations", 100);
		_logfile->list("Will run up to " + toString(_numNRIterations) + " Newton-Raphson iterations. (parameter 'numNRIterations)");
		_logfile->endIndent();
		_estimateExponential = true;
	}

	_pmdTables.initialize(&_bamFile._readGroups, _maxLengthForInference, _bamFile.maxReadLength(), _readGroupMap);
};

TPMDEstimator::~TPMDEstimator(){
	delete _readGroupMap;
};

void TPMDEstimator::_handleAlignment(){
	//check if it is forward or reverse strand!
	if(_alignment.isReverseStrand()){
		int d = 0;
		for(auto& b : _alignment){
			if(b.isAligned() && b.base != N){
				Base ref = _genoMap.flipBase(_alignment.referenceAtInternalPos(d));
				Base read = _genoMap.baseToFlippedBase[b.base];
				_pmdTables.addFromFivePrime(b.readGroupID, b.distFrom5Prime, ref, read);
				_pmdTables.addFromThreePrime(b.readGroupID, b.distFrom3Prime, ref, read);
			}
			++d;
		}
	} else {
		int d = 0;
		for(auto& b : _alignment){
			if(b.isAligned() && b.base != N){
				Base ref = _genoMap.toBase(_alignment.referenceAtInternalPos(d));
				_pmdTables.addFromFivePrime(b.readGroupID, b.distFrom5Prime, ref, b.base);
				_pmdTables.addFromThreePrime(b.readGroupID, b.distFrom3Prime, ref, b.base);
			}
			++d;
		}
	}
};

void TPMDEstimator::estimatePMD(){
	//traverse BAM
	_traverseBAMPassedQC();

	//print tables and data
	std::string filename = _outputName + "_PMD_Table.txt";
	_logfile->listFlush("Writing PMD table to '" + filename + "' ...");
	_pmdTables.writeTable(filename);
	_logfile->done();
	filename = _outputName + "_PMD_Table_counts.txt";
	_logfile->listFlush("Writing PMD table of counts to '" + filename + "' ...");
	_pmdTables.writeTableWithCounts(filename);
	_logfile->done();
	filename = _outputName + "_PMD_input_Empiric.txt";
	_logfile->listFlush("Writing PMD input file to '" + filename + "' ...");
	_pmdTables.writePMDFile(filename);
	_logfile->done();

	//estimate exponential model
	if(_estimateExponential){
		filename = _outputName + "_PMD_input_Exponential.txt";
		_logfile->listFlush("Estimating PMD exponential models and writing them to '" + filename + "' ...");
		_pmdTables.fitExponentialModel(_numNRIterations, _eps, filename, _logfile);
		_logfile->done();
	}
};

}; // end namespace
