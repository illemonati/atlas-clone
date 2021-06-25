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
	//make sure there is pmd
	if(!_genotypeLikelihoodCalculator.hasPMD()){
		throw "Can not estimate PMD: no PMD models provided! Use 'pmd' to provide PMD models.";
	}

	//make sure it has a reference
	_openReference(true);

	//prepare maps
	_readGroupMap = new BAM::TReadGroupMap(_bamFile.readGroups(), Parameters.getParameter<std::string>("poolReadGroups", false), _logfile);

	//parse estimation parameters
	_logfile->startIndent("Parameters for PMD Estimation:");
	_maxLengthForInference = Parameters.getParameterWithDefault<int>("length", 50);
	_logfile->list("Estimating PMD from the first " + coretools::str::toString(_maxLengthForInference) + " positions.");

	GenotypeLikelihoods::TPostMortemDamage& pmd = _genotypeLikelihoodCalculator.getPostMortemDamageModelsMutable();
	for(auto& r : _readGroupMap->readGroupsInUse()){
		pmd[r].parseEstimationParameters(_estimationParameters, Parameters, _logfile);
	}
	_logfile->endIndent();

	//create PMD tables
	_pmdTables.initialize(&_bamFile.readGroups(), _maxLengthForInference, _readGroupMap);
};

TPMDEstimator::~TPMDEstimator(){
	delete _readGroupMap;
};

void TPMDEstimator::_handleAlignment(){
	uint32_t d = 0;
	for(auto& b : _alignment){
		genometools::Base ref = _alignment.referenceAtInternalPos(d);
		_pmdTables.add(b, ref);
		++d;
	}
};

void TPMDEstimator::estimatePMD(){
	// 1) traverse BAM to assemble PMD Tables
	_traverseBAMPassedQC();

	// 2) print PMD tables
	_logfile->startIndent("Writing PMD tables:");

	//counts
	std::string filename = _outputName + "_PMD_counts.txt";
	_logfile->listFlush("Writing PMD counts to '" + filename + "' ...");
	_pmdTables.write(filename, false);
	_logfile->done();

	//normalized counts
	filename = _outputName + "_PMD_countsNormalized.txt";
	_logfile->listFlush("Writing PMD normalized counts to '" + filename + "' ...");
	_pmdTables.write(filename, true);
	_logfile->done();
	_logfile->endIndent();

	// 3) estimate models
	GenotypeLikelihoods::TPostMortemDamage& pmd = _genotypeLikelihoodCalculator.getPostMortemDamageModelsMutable();

	//estimate all models with data, i.e. only one model per pool
	for(auto& r : _readGroupMap->readGroupsInUse()){
		pmd[r].estimate(_pmdTables[r], _estimationParameters);
	}

	//writing PMD file
	filename = _outputName + "_ PMD.txt";
	pmd.writeToFile(_bamFile.readGroups(), *_readGroupMap, filename);
};

}; // end namespace
