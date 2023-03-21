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
TPMDEstimator::TPMDEstimator(): TGenome_parsed(), _readGroupMap(_bamFile.readGroups(), parameters().getParameter<std::string>("poolReadGroups", false)) {
	//make sure there is pmd
	GenotypeLikelihoods::TPostMortemDamage& pmd = _genotypeLikelihoodCalculator.postMortemDamageModels();
	if (_genotypeLikelihoodCalculator.hasPMD() && !parameters().parameterExists("reestimate")) {
		logfile().list("PMD model already exists, will reestimate it.");
	}
	if (!_genotypeLikelihoodCalculator.hasPMD()) {
		pmd.initialize(parameters().getParameterWithDefault("pmd", "doubleStrand:Empiric:Empiric"s), _bamFile.readGroups());
	}

	//make sure it has a reference
	_openReference(true);

	//parse estimation parameters
	logfile().startIndent("Parameters for PMD Estimation:");
	_maxLengthForInference = parameters().getParameterWithDefault<int>("length", 50);
	logfile().list("Estimating PMD from the first ", _maxLengthForInference, " positions.");

	//create PMD tables
	_pmdTables.initialize(&_bamFile.readGroups(), _maxLengthForInference, &_readGroupMap);
};

void TPMDEstimator::_handleAlignment() {
	for (size_t d = 0; d < _alignment.size(); ++d) {
		if (_alignment.isAlignedAtInternalPos(d)) {
			_pmdTables.add(_alignment[d], _alignment.referenceAtInternalPos(d));
		}
	}
};

void TPMDEstimator::run(){
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
	GenotypeLikelihoods::TPostMortemDamage& pmd = _genotypeLikelihoodCalculator.postMortemDamageModels();

	//estimate all models with data, i.e. only one model per pool
	pmd.estimate(_readGroupMap, _pmdTables);

	//writing PMD file
	filename = _outputName + "_PMD.txt";
	pmd.writeToFile(_bamFile.readGroups(), _readGroupMap, filename);
};

}; // end namespace
