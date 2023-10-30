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
using coretools::instances::parameters;

//----------------------------------------
// TPMDEstimator.h
//----------------------------------------
TPMDEstimator::TPMDEstimator()
	: TGenome_parsed(),
	  _readGroupMap(_bamFile.readGroups(), parameters().get("poolReadGroups", "")) {
	_pmd.initialize(parameters().get("pmdModel", "doubleStrand:Empiric:Empiric"),
					 _bamFile.readGroups());

	// make sure it has a reference
	_openReference(true);

	// create PMD tables
	_pmd.resize(_readGroupMap);
};

void TPMDEstimator::_handleAlignment() {
	for (size_t d = 0; d < _alignment.size(); ++d) {
		if (_alignment.isAlignedAtInternalPos(d)) {
			const auto data = _alignment[d];
			const auto from = _alignment.referenceAtInternalPos(d);
			_pmd.add(_readGroupMap, from, data);
		}
	}
};

void TPMDEstimator::run(){
	_traverseBAMPassedQC();
	_pmd.estimate(_readGroupMap);
	_pmd.writeToFile(_bamFile.readGroups(), _readGroupMap, _outputName);
};

}; // end namespace
