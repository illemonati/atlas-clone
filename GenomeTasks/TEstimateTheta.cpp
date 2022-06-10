/*
 * TEstimateTheta.cpp
 *
 *  Created on: Jun 4, 2020
 *      Author: phaentu
 */

#include "TEstimateTheta.h"

#include <algorithm>
#include <exception>
#include <map>
#include <memory>
#include <stddef.h>

#include "TBamFile.h"
#include "TFile.h"
#include "TGenomePosition.h"
#include "TGenotypeData.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TThetaEstimatorData.h"
#include "TTimer.h"
#include "stringFunctions.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::instances::randomGenerator;
using coretools::str::toString;

//-----------------------------------
// TEstimateTheta_base
//-----------------------------------
TEstimateTheta_base::TEstimateTheta_base()
	: TGenome_windows(), _thetaEstimator(){

						 };

void TEstimateTheta_base::_addSites(GenotypeLikelihoods::TWindow_base &window,
									GenotypeLikelihoods::TThetaEstimator &thetaEstimator) {
	logfile().listFlushTime("Calculating genotype likelihoods ...");
	for (auto &s : window) {
		_genoLik = _genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(s);
		thetaEstimator.add(s, _genoLik);
	}
	logfile().doneTime();
};

void TEstimateTheta_base::_addSites() { _addSites(_window, _thetaEstimator); };

//-----------------------------------
// TEstimateTheta
//-----------------------------------
TEstimateTheta::TEstimateTheta() : TEstimateTheta_base() {
	// open output file
	std::string filename = _outputName + "_thetaPerWindow.txt.gz";
	_thetaOut.open(&_thetaEstimator, filename);

	// print all windows?
	if (parameters().parameterExists("printAll")) {
		_printAll = true;
		logfile().list(
			"Will print all windows, also those for which no estimation was possible. (parameter 'printAll')");
	} else {
		_printAll = false;
		logfile().list("Will only print windows for which estimation was possible. (use 'printAll' to print all)");
	}
};

bool TEstimateTheta::_estimateTheta() {
	_addSites();

	// estimate Theta
	return _thetaEstimator.estimateTheta();
};

void TEstimateTheta::_handleWindow() {
	logfile().startIndent("Estimating Theta:");
	if (_estimateTheta() || _printAll) { _thetaOut.write(_window); }
	_thetaEstimator.clear();
	logfile().endIndent();
};

void TEstimateTheta::estimateTheta() { _traverseBAMWindows(); };

//-----------------------------------
// TEstimateThetaGenomeWide
//-----------------------------------
TEstimateThetaGenomeWide::TEstimateThetaGenomeWide() : TEstimateTheta_base() {
	if (_considerRegions) {
		logfile().list("Estimating theta at specific sites. (parameter 'regions')");
	} else {
		logfile().list("Estimating theta genome-wide. (use 'regions' to limit)");
	}

	// bootstraps
	_numBootstraps = parameters().getParameterWithDefault<int>("bootstraps", 0);
	if (_numBootstraps > 0) {
		logfile().list("Will estimate theta fpr ", _numBootstraps, " bootstrap replicates. (parameter 'bootstraps')");
	} else {
		logfile().list("Will not conduct any bootstrap replicates. (use 'bootstraps' to request)");
	}

	if (parameters().parameterExists("onlyBootstrap")) {
		_onlyBootstraps = true;
		logfile().list("Will only ...");
	} else {
		_onlyBootstraps = true;
		logfile().list("Will only ...");
	}
};

void TEstimateThetaGenomeWide::_handleWindow() {
	// add sites to data structure
	try {
		_addSites();
	} catch (...) {
		throw "Failed to allocate sufficient memory to store the data for so many sites. Consider limiting the "
			  "analysis (e.g. parameters 'regions'or 'chr') or limiting to sites with depth >=2 (parameter 'minDepth', "
			  "recommended).";
	}
};

void TEstimateThetaGenomeWide::_bootstrapThetaEstimation() {
	logfile().startIndent("Generating " + toString(_numBootstraps) + " bootstrap estimates of theta:");

	// measure runtime
	coretools::TTimer timer;

	// loop over bootstraps
	for (uint32_t s = 0; s < _numBootstraps; ++s) {
		logfile().startIndent("Bootstrap " + toString(s + 1) + " of " + toString(_numBootstraps) + ":");

		// run bootstrap
		_thetaEstimator.bootstrapTheta();
		_thetaOut.write("Bootstrap_" + toString(s + 1), "-", "-");

		logfile().endIndent();
	}

	// finish
	logfile().list("Total computation time for theta bootstrapping was ", timer.minutes());
	logfile().endIndent();
};

void TEstimateThetaGenomeWide::estimateThetaGenomeWide() {
	// read data
	_traverseBAMWindows();

	if (!_onlyBootstraps) {
		logfile().startIndent("Estimate theta based on a total of " + toString(_thetaEstimator.sizeWithData()) +
							  " sites:");
		_thetaEstimator.estimateTheta();

		// write estimates
		std::string filename = _outputName + "_thetaGenomeWide.txt.gz";
		_thetaOut.open(&_thetaEstimator, filename);
		logfile().list("Will write theta estimates to file '" + filename + "'.");
		if (_considerRegions) {
			_thetaOut.write("regions", "-", "-");
		} else {
			_thetaOut.write("genome-wide", "-", "-");
		}
		logfile().endIndent();
	}

	// bootstrap
	if (_numBootstraps > 0) { _bootstrapThetaEstimation(); }
};

//-----------------------------------
// TEstimateThetaLLSurface
//-----------------------------------
TEstimateThetaLLSurface::TEstimateThetaLLSurface() : TEstimateTheta_base() {
	_steps = parameters().getParameterWithDefault<int>("steps", 100);
	logfile().list("Will calculate the LL-surface at ", _steps, " steps. (parameter 'steps')");
	if (_steps < 2) { throw "Th enumber of steps must be >= 2!"; }
};

void TEstimateThetaLLSurface::_handleWindow() {
	logfile().startIndent("Calculating likelihood surface for Theta:");

	// adding sites to estimator
	_addSites();

	// open file
	std::string filename =
		_outputName + _window.chrName() + "_" + toString(_window.from().position()) + "_LLsurface.txt";
	logfile().listFlushTime("Writing LL surface to file '" + filename + "' ...");
	coretools::TOutputFile out(filename);

	_thetaEstimator.calcLikelihoodSurface(out, _steps);

	_thetaEstimator.clear();
	logfile().doneTime();

	// finish
	logfile().endIndent();
};

void TEstimateThetaLLSurface::estimateThetaLLSurface() { _traverseBAMWindows(); };

//-----------------------------------
// TEstimateThetaDownsamplingQC
//-----------------------------------
TEstimateThetaDownsamplingQC::TEstimateThetaDownsamplingQC() : TEstimateTheta_base() {

	// read downsampling rates
	if (parameters().parameterExists("prob")) {
		parameters().fillParameterIntoContainer("prob", downSampleProbVector, ',');
	} else if (parameters().parameterExists("depth")) {
		std::vector<double> depths;
		parameters().fillParameterIntoContainer("depth", depths, ',');
		double averageDepth = parameters().getParameter<double>("averageDepth");
		for (auto &it : depths) {
			if (averageDepth >= it) {
				downSampleProbVector.push_back(it / averageDepth);
			} else {
				throw "Average Depth must be equal or bigger than provided lists of depths";
			}
		}
	} else {
		// parse downsampling parameters
		parameters().fillParameterIntoContainer("prob", downSampleProbVector, ',', "1.0,0.5,0.2,0.1,0.05,0.02,0.01");
	}
	// report probabilities
	logfile().list("Will estimate theta after downsampling reads with probabilities " +
				   coretools::str::concatenateString(downSampleProbVector, ", ") + ". (parameter 'prob')");

	// check if full data is to be used (i.e. if prob = 1.0 is specified)
	_printFullData = false;
	if (*downSampleProbVector.begin() == 1.0) {
		_printFullData = true;
		downSampleProbVector.erase(downSampleProbVector.begin());
	}

	// create windows and estimators for downsamples
	if (downSampleProbVector.size() > 0) {
		for (size_t i = 0; i < downSampleProbVector.size(); ++i) { estimators.emplace_back(_thetaEstimator); }
	}

	// open output
	std::string filename = _outputName + "_thetaQC.txt.gz";
	logfile().list("Will write theta estimates to file '" + filename + "'.");
	if (_printFullData) { _thetaOut.addEstimator(&_thetaEstimator, "p1.0_"); }
	for (size_t i = 0; i < downSampleProbVector.size(); ++i) {
		// assemble prefix without lagging zeros
		std::string prefix = toString(downSampleProbVector[i]);
		int pos            = prefix.size() - 1;
		while (prefix[pos] == '0') pos--;
		prefix.erase(pos + 1, prefix.size() - pos - 1);
		prefix = "p" + prefix + "_";

		// now add estimator to output file
		_thetaOut.addEstimator(&estimators[i], prefix);
	}
	_thetaOut.open(filename);
};

void TEstimateThetaDownsamplingQC::_handleWindow() {
	// estimate on full data
	if (_printFullData) {
		logfile().startIndent("Estimating Theta on full data:");
		_addSites();
		_thetaEstimator.estimateTheta();
		logfile().endIndent();
	}

	for (size_t i = 0; i < downSampleProbVector.size(); ++i) {
		coretools::Probability &p = downSampleProbVector[i];
		logfile().startIndent("Estimating Theta on downsampled data (p = " + toString(p) + "):");

		// downsample
		logfile().listFlush("Downsampling reads with p = " + toString(p) + " ...");
		if (_subset)
			destination.downsampleFromOther(_window, *_subset, _readUpToDepth, p, &randomGenerator());
		else
			destination.downsampleFromOther(_window, _readUpToDepth, p, &randomGenerator());
		logfile().done();

		// apply filters
		_applyWindowFilters(destination);

		// and estimate
		_addSites(destination, estimators[i]);
		estimators[i].estimateTheta();
		logfile().endIndent();
	}

	// write output
	_thetaOut.write(_window);

	// clear
	_thetaEstimator.clear();
	for (auto &e : estimators) { e.clear(); }
};

void TEstimateThetaDownsamplingQC::runQC() { _traverseBAMWindows(); };

//-----------------------------------
// TEstimateThetaRatio
//-----------------------------------
TEstimateThetaRatio::TEstimateThetaRatio() : TGenome_windows(), _thetaEstimatorRatio() {
	// read the two regions to be used
	logfile().startIndent("Reading regions:");
	_initializeRegion(_region1, '1');
	_initializeRegion(_region1, '2');
};

void TEstimateThetaRatio::_initializeRegion(genometools::TBed &region, const char num) {
	logfile().startIndent((std::string) "Region " + num + ":");
	std::string regionsFile = parameters().getParameter<std::string>("regions" + std::to_string(num));
	logfile().list((std::string) "Reading regions " + num + " from file '" + regionsFile + " (parameter 'region" + num +
				   "') ...");
	region.add(regionsFile, _bamFile.chromosomes());
	logfile().done();
	logfile().conclude("Read " + toString(region.size()) + " sites on " + toString(region.numChromosomesWithWindows()) +
					   " chromosomes.");
};

void TEstimateThetaRatio::_addSites(GenotypeLikelihoods::TThetaEstimatorData &data, genometools::TBed &region) {
	auto it = region.lower_bound(_window);

	while (it != region.end() && _window.overlaps(*it)) {
		for (genometools::TGenomePosition s = std::max(it->from(), _window.from()); s < it->to() && s < _window.to();
			 ++s) {
			GenotypeLikelihoods::TGenotypeLikelihoods genoLik;
			genoLik = _genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(_window[s - _window.from()]);
			data.add(_window[s - _window.from()], genoLik);
		}
		++it;
	}
};

void TEstimateThetaRatio::_handleWindow() {
	// adding sites to estimator
	logfile().listFlushTime("Calculating genotype likelihoods ...");
	try {
		_addSites(*_thetaEstimatorRatio.pointerToDataContainer(), _region1);
		_addSites(*_thetaEstimatorRatio.pointerToDataContainer2(), _region2);
	} catch (...) {
		throw "Failed to allocate sufficient memory to store the data for so many sites. Consider selecting fewer "
			  "regions or limiting to sites with a minimal depth (>=2 recommended).";
	}
	logfile().doneTime();
};

void TEstimateThetaRatio::estimateThetaRation() {
	_traverseBAMWindows();
	_thetaEstimatorRatio.estimateRatio(_outputName);
};

}; // namespace GenomeTasks
