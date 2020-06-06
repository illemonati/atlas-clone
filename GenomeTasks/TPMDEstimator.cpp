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
TPMDEstimator::TPMDEstimator(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_parsed(Params, Logfile, RandomGenerator){
	//prepare maps
	_readGroupMap = new BAM::TReadGroupMap(&_bamFile.readGroups, Params.getParameterString("poolReadGroups", false), _logfile);

	//prepare PMD table
	_maxLengthForInference = Params.getParameterIntWithDefault("length", 50);
	_logfile->list("Estimating PMD at the first " + toString(_maxLengthForInference) + " positions.");

	if(Params.parameterExists("onlyEmpiric")){
		_logfile->list("Not fitting exponential models. (parameter 'onlyEmpiric')");
		_estimateExponential = false;
	} else {

		_logfile->startIndent("Will estimate exponential PMD models: (use 'onlyEmpiric' to turn off)");
		_eps = Params.getParameterDoubleWithDefault("eps", 0.001);
		_logfile->list("Will conider the Newton-Raphson algorithm to have converged if the likelihood difference < " + toString(_eps) + ". (parameter 'eps')");
		_numNRIterations = Params.getParameterIntWithDefault("numNRIterations", 100);
		_logfile->list("Will run up to " + toString(_numNRIterations) + " Newton-Raphson iterations. (parameter 'numNRIterations)");
		_logfile->endIndent();
		_estimateExponential = true;
	}

	_pmdTables.initialize(&_bamFile.readGroups, _maxLengthForInference, _bamFile.maxReadLength(), _readGroupMap);
};

TPMDEstimator::~TPMDEstimator(){
	delete _readGroupMap;
};

void TPMDEstimator::_handleAlignment(){
	_alignment.addToPMDTables(_pmdTables, _genoMap);
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
