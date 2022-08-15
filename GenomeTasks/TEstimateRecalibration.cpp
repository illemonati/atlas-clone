/*
 * TestimateRecalibration.cpp
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#include "TEstimateRecalibration.h"

#include <exception>
#include <memory>
#include <stdint.h>
#include <vector>

#include "TBamFile.h"
#include "TGenomePosition.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TLog.h"
#include "TParameters.h"
#include "TReadGroups.h"
#include "SequencingError/TRecalEstimator.h"
#include "TSite.h"
#include "TSiteSubset.h"
#include "TWindow.h"

namespace GenomeTasks {
using coretools::instances::logfile;
using coretools::instances::parameters;

//-----------------------------------------------------------
// TEstimateRecalibration_base
//-----------------------------------------------------------
TEstimateRecalibration::TEstimateRecalibration() : TGenome_windows() {
	if (_genotypeLikelihoodCalculator.recalibrationChangesQualities() && !parameters().parameterExists("rerecalibrate"))
		throw "Can not estimate recalibration: quality scores are already recalibrated while reading! (Use argument "
			  "'rerecalibrate' to overwrite this error)";

	// limit to sites with known alleles?
	if (parameters().parameterExists("alleles")) {
		logfile().startIndent("Will limit analysis to sites with known genotypes (parameter 'genotypes'):");
		_openSiteSubset("alleles");
		logfile().endIndent();
	} else {
		logfile().list(
			"Will use all sites without prior knowledge on alleles. (use 'genotypes' to provide known genotypes)");
	}

	// initialize maps
	if (parameters().parameterExists("poolReadGroups")) {
		std::string poolReadGroupsFile = parameters().getParameter<std::string>("poolReadGroups");
		logfile().startIndent("Will pool read groups (parameter 'poolReadGroups'):");
		_readGroupMap = std::make_unique<BAM::TReadGroupMap>(_bamFile.readGroups(), poolReadGroupsFile, &logfile());
		logfile().endIndent();
	} else {
		logfile().list("Will estimate recalibration parameters for each read group. (use 'poolReadGroups' to pool)");
		_readGroupMap = std::make_unique<BAM::TReadGroupMap>(_bamFile.readGroups());
	}

	_onlyLL = parameters().parameterExists("onlyLL");

	// initialize recal estimator
	recalObjectEM = std::make_unique<GenotypeLikelihoods::SequencingError::TRecalibrationEMEstimator>(
		&(_bamFile.readGroups()), _readGroupMap.get());
};

void TEstimateRecalibration::_handleWindow() {
	// add sites to recal estimator
	if (_subset) {
		std::set<GenotypeLikelihoods::TSiteSubsetSite> thesePositions = _subset->getPositionInWindow(_window);
		for (auto &it : thesePositions) {
			uint32_t internalPos = it - _window.from();
			if (!_window[internalPos].empty() && it.ref() == it.alt()) { recalObjectEM->addSite(_window[internalPos]); }
		}
	} else {
		for (auto &s : _window) {
			if (!s.empty()) { recalObjectEM->addSite(s); }
		}
	}
};

void TEstimateRecalibration::estimateRecalibration() {
	// read data
	_traverseBAMWindows();

	if (_onlyLL) {
		recalObjectEM->calcLL(_genotypeLikelihoodCalculator.getSequencingErrorModelsMutable(),
							  _genotypeLikelihoodCalculator.getPostMortemDamageModelsMutable());

	} else {
		// estimate recal parameters
		recalObjectEM->performEstimation(_outputName, _genotypeLikelihoodCalculator.getSequencingErrorModelsMutable(),
										 _genotypeLikelihoodCalculator.getPostMortemDamageModelsMutable());
	}
};

}; // namespace GenomeTasks
