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
	_readGroupMap = new BAM::TReadGroupMap(&_bamFile.readGroupsMutable(), Parameters.getParameterString("poolReadGroups", false), _logfile);

	//prepare PMD table
	_maxLengthForInference = Parameters.getParameterIntWithDefault("length", 50);
	_logfile->list("Estimating PMD from the first " + toString(_maxLengthForInference) + " positions.");
	_pmdTables.initialize(&_bamFile.readGroupsMutable(), _maxLengthForInference, _bamFile.maxReadLength(), _readGroupMap);

	//create PMD objects
	_pmd.initialize(Parameters, _bamFile.readGroupsMutable(), Logfile);
};

TPMDEstimator::~TPMDEstimator(){
	delete _readGroupMap;
};

void TPMDEstimator::_handleAlignment(){
	uint32_t d = 0;
	for(auto& b : _alignment){
		Base ref = _genoMap.flipBase(_alignment.referenceAtInternalPos(d));
		_pmdTables.add(b, ref);
		++d;
	}
};

void TPMDEstimator::estimatePMD(){
	//check if PMD must be estimated
	if(!_pmd.hasPMD()){
		_logfile->list("Data has no PMD. Aborting estimation.");
	} else {
		//estimate PMD
		//traverse BAM
		_traverseBAMPassedQC();

		//print PMD tables
		_logfile->startIndent("Writing PMD tables:");

		//counts
		filename = _outputName + "_PMD_counts.txt";
		_logfile->listFlush("Writing PMD counts to '" + filename + "' ...");
		_pmdTables.write(filename, false);
		_logfile->done();

		//normalized counts
		std::string filename = _outputName + "_PMD_countsNormalized.txt";
		_logfile->listFlush("Writing PMD normalized counts to '" + filename + "' ...");
		_pmdTables.writeTable(filename);
		_logfile->done();
		_logfile->endIndent();

		//estimate models
		_logfile->startIndent("Estimating PMD functions:");


			_logfile->listFlush("Estimating PMD exponential models and writing them to '" + filename + "' ...");
			_pmdTables.fitExponentialModel(_numNRIterations, _eps, filename, _logfile);
			_logfile->done();

			_logfile->endIndent();
		}
	}

	//writing PMD file
	std::string filename = _outputName + "_ PMD.txt";
	_pmd.writeToFile(_bamFile.readGroupsMutable(), filename);
};

}; // end namespace
