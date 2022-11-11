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

#include "genometools/GenotypeTypes.h"
#include "TAlignment.h"
#include "TBamFile.h"
#include "TGenome.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "TReadGroups.h"

namespace GenomeTasks{
using namespace std::literals;
using coretools::instances::logfile;
using coretools::instances::parameters;

//----------------------------------------
// TPMDEstimator.h
//----------------------------------------
TPMDEstimator::TPMDEstimator(): TGenome_parsed() {
	//make sure there is pmd
	GenotypeLikelihoods::TPostMortemDamage& pmd = _genotypeLikelihoodCalculator.getPostMortemDamageModels();
	if (_genotypeLikelihoodCalculator.hasPMD() && !parameters().parameterExists("reestimate")) {
		throw "PMD model already estimated! (Use argument 'reestimate' to overwrite this error)";
	}
	if (!_genotypeLikelihoodCalculator.hasPMD()) {
		pmd.initialize(parameters().getParameterWithDefault("pmdModels", "doubleStrand:Empiric:Empiric"s), _bamFile.readGroups());
	}

	//make sure it has a reference
	_openReference(true);

	//prepare maps
	_readGroupMap = std::make_unique<BAM::TReadGroupMap>(_bamFile.readGroups(), parameters().getParameter<std::string>("poolReadGroups", false), &logfile());

	//parse estimation parameters
	logfile().startIndent("Parameters for PMD Estimation:");
	_maxLengthForInference = parameters().getParameterWithDefault<int>("length", 50);
	logfile().list("Estimating PMD from the first ", _maxLengthForInference, " positions.");

	for(auto& r : _readGroupMap->readGroupsInUse()){
		pmd[r].parseEstimationParameters(_estimationParameters);
	}
	logfile().endIndent();

	//create PMD tables
	_pmdTables.initialize(&_bamFile.readGroups(), _maxLengthForInference, _readGroupMap.get());
};

void TPMDEstimator::_handleAlignment() {
	for (size_t d = 0; d < _alignment.size(); ++d) {
		if (_alignment.isAlignedAtInternalPos(d)) {
			_pmdTables.add(_alignment[d], _alignment.referenceAtInternalPos(d));
		}
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
	GenotypeLikelihoods::TPostMortemDamage& pmd = _genotypeLikelihoodCalculator.getPostMortemDamageModels();

	//estimate all models with data, i.e. only one model per pool
	for(auto& r : _readGroupMap->readGroupsInUse()){
		pmd[r].estimate(_pmdTables[r], _estimationParameters);
	}

	//writing PMD file
	filename = _outputName + "_PMD.txt";
	pmd.writeToFile(_bamFile.readGroups(), *_readGroupMap, filename);
};

}; // end namespace
