/*
 * TestimateRecalibration.cpp
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#include "TEstimateErrors.h"
#include "PMD/TModels.h"


namespace GenomeTasks {
using coretools::instances::logfile;
using coretools::instances::parameters;

namespace impl {
BAM::TReadGroupMap makeRGMap(const BAM::TReadGroups &ReadGroups) {
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
// TEstimateErrors_base
//-----------------------------------------------------------
TEstimateErrors::TEstimateErrors()
    : TGenome_windows(), _errors(_bamFile.readGroups()) {
	_openReference(true);

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

void TEstimateErrors::_handleWindow() {
	// add sites to recal estimator
	if (_subsetMonomorphic) {
		auto thesePositions = _subsetMonomorphic->getPositionInWindow(_window);
		for (auto &it : thesePositions) {			
			const uint32_t internalPos = it - _window.from();
			_errors.addSite(_window[internalPos]);
		}
	} else {
		for (auto &s : _window) {
			_errors.addSite(s);
		}
	}
};

void TEstimateErrors::run() {
	// read data
	_traverseBAMWindows();

	if (_onlyLL)
		_errors.calcLL();
	else
		_errors.estimate(_outputName);
};

}; // namespace GenomeTasks
