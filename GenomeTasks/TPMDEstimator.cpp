/*
 * TPMDEstimator.cpp
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#include "TPMDEstimator.h"

#include <exception>
#include <memory>
#include <vector>

#include "GenotypeTypes.h"
#include "TAlignment.h"
#include "TBamFile.h"
#include "TGenome.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TLog.h"
#include "TParameters.h"
#include "TReadGroups.h"

namespace GenomeTasks{
using namespace std::literals;
using coretools::instances::logfile;
using coretools::instances::parameters;

//----------------------------------------
// TPMDEstimator.h
//----------------------------------------
TPMDEstimator::TPMDEstimator(): TGenome_parsed(){
	//make sure there is pmd
	if(!_genotypeLikelihoodCalculator.hasPMD()){
		throw "Can not estimate PMD: no PMD models provided! Use 'pmd' to provide PMD models.";
	}

	//make sure it has a reference
	_openReference(true);

	//prepare maps
	_readGroupMap = std::make_unique<BAM::TReadGroupMap>(_bamFile.readGroups(), parameters().getParameter<std::string>("poolReadGroups", false), &logfile());

	//parse estimation parameters
	logfile().startIndent("Parameters for PMD Estimation:");
	_maxLengthForInference = parameters().getParameterWithDefault<int>("length", 50);
	logfile().list("Estimating PMD from the first ", _maxLengthForInference, " positions.");

	GenotypeLikelihoods::TPostMortemDamage& pmd = _genotypeLikelihoodCalculator.getPostMortemDamageModelsMutable();
	for(auto& r : _readGroupMap->readGroupsInUse()){
		pmd[r].parseEstimationParameters(_estimationParameters);
	}
	logfile().endIndent();

	//create PMD tables
	_pmdTables.initialize(&_bamFile.readGroups(), _maxLengthForInference, _readGroupMap.get());
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
	logfile().startIndent("Writing PMD tables:");

	//counts
	std::string filename = _outputName + "_PMD_counts.txt";
	logfile().listFlush("Writing PMD counts to '" + filename + "' ...");
	_pmdTables.write(filename, false);
	logfile().done();

	//normalized counts
	filename = _outputName + "_PMD_countsNormalized.txt";
	logfile().listFlush("Writing PMD normalized counts to '" + filename + "' ...");
	_pmdTables.write(filename, true);
	logfile().done();
	logfile().endIndent();

	// 3) estimate models
	using coretools::instances::parameters;
	GenotypeLikelihoods::TPostMortemDamage& pmd = _genotypeLikelihoodCalculator.getPostMortemDamageModelsMutable();
	if(pmd.hasPMD()) {
		if (!parameters().parameterExists("reestimate")) throw "not happy";
	} else {
		std::vector<uint16_t> tmp;
		pmd.initialize(parameters().getParameterWithDefault("pmdModels", "doubleStrand:Empirical:Empirical"s), _bamFile.readGroups(), tmp);
	}

	//estimate all models with data, i.e. only one model per pool
	for(auto& r : _readGroupMap->readGroupsInUse()){
		pmd[r].estimate(_pmdTables[r], _estimationParameters);
	}

	//writing PMD file
	filename = _outputName + "_PMD.txt";
	pmd.writeToFile(_bamFile.readGroups(), *_readGroupMap, filename);
};

}; // end namespace
