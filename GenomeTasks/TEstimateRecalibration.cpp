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
#include "genometools/GenomePositions/TGenomePosition.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "TReadGroups.h"
#include "SequencingError/TRecalEstimator.h"
#include "TSite.h"
#include "TSiteSubset.h"
#include "TWindow.h"

namespace GenomeTasks {
using coretools::instances::logfile;
using coretools::instances::parameters;

namespace impl {
BAM::TReadGroupMap makeReadGroupMap(const BAM::TReadGroups &ReadGroups) {
	if (parameters().parameterExists("poolReadGroups")) {
		std::string poolReadGroupsFile = parameters().getParameter<std::string>("poolReadGroups");
		logfile().startIndent("Will pool read groups (parameter 'poolReadGroups'):");
		return {ReadGroups, poolReadGroupsFile};
		logfile().endIndent();
	} else {
		logfile().list("Will estimate recalibration parameters for each read group. (use 'poolReadGroups' to pool)");
		return {ReadGroups};
	}
}
} // namespace impl

//-----------------------------------------------------------
// TEstimateRecalibration_base
//-----------------------------------------------------------
	TEstimateRecalibration::TEstimateRecalibration() : TGenome_windows(),
													   _readGroupMap(impl::makeReadGroupMap(_bamFile.readGroups())),
													   recal(&_bamFile.readGroups(), &_readGroupMap) {
		if (_genotypeLikelihoodCalculator.recalibrationChangesQualities() && !(
				parameters().parameterExists("rerecalibrate")
				|| parameters().getParameter("task") == "recal"))
			UERROR("Can not estimate recalibration: quality scores are already recalibrated while reading! (Use argument "
				   "'rerecalibrate' to overwrite this error)");

	// limit to sites with known alleles?
	if (parameters().parameterExists("alleles")) {
		logfile().startIndent("Will limit analysis to sites with known genotypes (parameter 'genotypes'):");
		_openSiteSubset("alleles");
		logfile().endIndent();
	} else {
		logfile().list(
			"Will use all sites without prior knowledge on alleles. (use 'genotypes' to provide known genotypes)");
	}

	_onlyLL = parameters().parameterExists("onlyLL");
};

void TEstimateRecalibration::_handleWindow() {
	// add sites to recal estimator
	if (_subset) {
		auto thesePositions = _subset->getPositionInWindow(_window);
		for (auto &it : thesePositions) {
			if (it.ref() == it.alt()) {
				const uint32_t internalPos = it - _window.from();
				recal.addSite(_window[internalPos]);
			}
		}
	} else {
		for (auto &s : _window) {
			recal.addSite(s);
		}
	}
};

void TEstimateRecalibration::estimateRecalibration() {
	// read data
	auto forgottenModels = _genotypeLikelihoodCalculator.sequencingErrorModels().forget();
	_traverseBAMWindows();
	_genotypeLikelihoodCalculator.sequencingErrorModels().remember(forgottenModels);

	if (_onlyLL) {
		recal.calcLL(_genotypeLikelihoodCalculator.sequencingErrorModels(),
							  _genotypeLikelihoodCalculator.postMortemDamageModels());

	} else {
		// estimate recal parameters
		recal.performEstimation(_outputName, _genotypeLikelihoodCalculator.sequencingErrorModels(),
										 _genotypeLikelihoodCalculator.postMortemDamageModels());
	}
};

}; // namespace GenomeTasks
