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

#include "coretools/Files/TOutputFile.h"
#include "genometools/GenotypeTypes.h"
#include "TAlignment.h"
#include "TBamFile.h"
#include "TGenome.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "TReadGroups.h"

namespace GenomeTasks{
using coretools::instances::logfile;
using coretools::instances::parameters;

//----------------------------------------
// TPMDEstimator.h
//----------------------------------------
	TPMDEstimator::TPMDEstimator(): TGenome_parsed(), _readGroupMap(_bamFile.readGroups(), parameters().getParameter<std::string>("poolReadGroups", false)), _pmd(&_genotypeLikelihoodCalculator.postMortemDamageModels()) {
	//make sure there is pmd
	if (_genotypeLikelihoodCalculator.hasPMD()) {
		logfile().list("PMD model already exists, will reestimate it.");
	}
	if (!_genotypeLikelihoodCalculator.hasPMD()) {
		_pmd->initialize(parameters().getParameterWithDefault("pmd", "doubleStrand:Empiric:Empiric"), _bamFile.readGroups());
	}

	//make sure it has a reference
	_openReference(true);

	//create PMD tables
	_pmd->resize(_readGroupMap);
};

void TPMDEstimator::_handleAlignment() {
	for (size_t d = 0; d < _alignment.size(); ++d) {
		if (_alignment.isAlignedAtInternalPos(d)) {
			const auto data     = _alignment[d];
			const auto from     = _alignment.referenceAtInternalPos(d);
			_pmd->add(_readGroupMap, from, data);
		}
	}
};

void TPMDEstimator::run(){
	_traverseBAMPassedQC();

	GenotypeLikelihoods::TPostMortemDamage& pmd = _genotypeLikelihoodCalculator.postMortemDamageModels();
	pmd.estimate(_readGroupMap);

	pmd.writeToFile(_bamFile.readGroups(), _readGroupMap, _outputName);
};

}; // end namespace
