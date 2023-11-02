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
#include "coretools/Files/TFile.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "TGenotypeData.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TThetaEstimatorData.h"
#include "coretools/TTimer.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::str::toString;


//-----------------------------------
// TEstimateThetaLLSurface
//-----------------------------------
TEstimateThetaLLSurface::TEstimateThetaLLSurface() : TGenome_windows() {
	_steps = parameters().getParameterWithDefault<int>("steps", 100);
	logfile().list("Will calculate the LL-surface at ", _steps, " steps. (parameter 'steps')");
	if (_steps < 2) { UERROR("Th enumber of steps must be >= 2!"); }
};

void TEstimateThetaLLSurface::_handleWindow() {
	logfile().startIndent("Calculating likelihood surface for Theta:");

	// adding sites to estimator
	for (auto &s : _window) {
		_genoLik = _genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(s);
		_thetaEstimator.add(s, _genoLik);
	}

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

void TEstimateThetaLLSurface::run() { _traverseBAMWindows(); };

//-----------------------------------
// TEstimateThetaDownsamplingQC
//-----------------------------------


void TEstimateTheta::_addSites(GenotypeLikelihoods::TWindow &window,
									GenotypeLikelihoods::TThetaEstimator &thetaEstimator) {
	logfile().listFlushTime("Calculating genotype likelihoods ...");
	for (auto &s : window) {
		_genoLik = _genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(s);
		thetaEstimator.add(s, _genoLik);
	}
	logfile().doneTime();
};

void TEstimateTheta::_addSites() { _addSites(_window, _thetaEstimator); };

TEstimateTheta::TEstimateTheta() : TGenome_windows() {
	if (parameters().parameterExists("genomeWide")) {
		_genomeWide = true;
		logfile().list("Will estimating heterozygosity (theta) genome-wide.");

		if (_considerRegions) {
			logfile().list("Estimating theta at specific sites. (parameter 'regions')");
		} else {
			logfile().list("Estimating theta genome-wide. (use 'regions' to limit)");
		}

		// bootstraps
		_numBootstraps = parameters().getParameterWithDefault<int>("bootstraps", 0);
		if (_numBootstraps > 0) {
			logfile().list("Will estimate theta fpr ", _numBootstraps,
						   " bootstrap replicates. (parameter 'bootstraps')");
		} else {
			logfile().list("Will not conduct any bootstrap replicates. (use 'bootstraps' to request)");
		}

		if (parameters().parameterExists("onlyBootstrap")) {
			_onlyBootstraps = true;
			logfile().list("Will only bootstrap");
		} else {
			_onlyBootstraps = false;
		}

	} else {
		_genomeWide = false;
		logfile().list("Will estimating heterozygosity (theta) per Window.");
	}

	// read downsampling rates

	if (parameters().parameterExists("prob")) {
		parameters().fillParameterIntoContainerRepeatIndexes("prob", downSampleProbVector, ',');
	} else if (parameters().parameterExists("depth")) {
		std::vector<double> depths;
		parameters().fillParameterIntoContainerRepeatIndexes("depth", depths, ',');
		double averageDepth = parameters().getParameter<double>("averageDepth");
		for (auto &it : depths) {
			if (averageDepth >= it) {
				downSampleProbVector.push_back(it / averageDepth);
			} else {
				UERROR("Average Depth must be equal or bigger than provided lists of depths");
			}
		}
	} else {
		downSampleProbVector.emplace_back(1.);
	}

	if (downSampleProbVector.empty()) UERROR("You need to specify at least one probability!");

	// check if full data is to be used (i.e. if prob = 1.0 is specified)
	_printFullData = false;
	if (downSampleProbVector.front() == 1.0) {
		_printFullData = true;
		downSampleProbVector.erase(downSampleProbVector.begin());
		logfile().list("Will estimate theta on full data.");
	}

	// create windows and estimators for downsamples
	if (downSampleProbVector.size() > 0) {
		logfile().list("Will estimate theta after downsampling reads with probabilities " +
					   coretools::str::concatenateString(downSampleProbVector, ", ") + ". (parameter 'prob')");
		for (size_t i = 0; i < downSampleProbVector.size(); ++i) { estimators.emplace_back(_thetaEstimator); }
	}

	// open output
	std::string filename = _outputName + "_theta.txt.gz";
	if (_printFullData) {
		const std::string prefix = downSampleProbVector.empty()? "" : "p1.0_";
		_thetaOut.addEstimator(&_thetaEstimator, prefix);
	}
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

	// print all windows?
	if (parameters().parameterExists("printAll")) {
		_printAll = true;
		logfile().list(
			"Will print all windows, also those for which no estimation was possible. (parameter 'printAll')");
	} else {
		_printAll = false;
		logfile().list("Will only print windows for which estimation was possible. (use 'printAll' to print all)");
	}
	logfile().endIndent();
};

void TEstimateTheta::_handleWindow() {
	// estimate on full data
	bool pass = false;
	if (_printFullData) {
		logfile().startIndent("Using full data:");

		_addSites();
		if (!_genomeWide) {
			logfile().startIndent("Estimating Theta:");

			pass |= _thetaEstimator.estimateTheta();
			logfile().endIndent();
		}
		logfile().endIndent();
	}

	static GenotypeLikelihoods::TWindow destination;

	for (size_t i = 0; i < downSampleProbVector.size(); ++i) {
		coretools::Probability &p = downSampleProbVector[i];
		logfile().startIndent("Using downsampled data (p = ", p, "):");

		logfile().listFlush("Downsampling reads ...");		
		destination.downsampleFromOther(_window, _readUpToDepth, p);
		logfile().done();

		_applyWindowFilters(destination);
		_addSites(destination, estimators[i]);

		if (!_genomeWide) {
			logfile().startIndent("Estimating Theta:");
			pass |= estimators[i].estimateTheta();
			logfile().endIndent();
		}
		logfile().endIndent();
	}

	// write output & clear
	if (!_genomeWide) {
		if (pass || _printAll) _thetaOut.write(_window);

		_thetaEstimator.clear();
		for (auto &e : estimators) { e.clear(); }
	}
};

void TEstimateTheta::_bootstrapThetaEstimation() {
	logfile().startIndent("Generating " + toString(_numBootstraps) + " bootstrap estimates of theta:");

	// measure runtime
	coretools::TTimer timer;

	// loop over bootstraps
	for (size_t s = 0; s < _numBootstraps; ++s) {
		logfile().startIndent("Bootstrap " + toString(s + 1) + " of " + toString(_numBootstraps) + ":");

		// run bootstrap
		_thetaEstimator.bootstrapTheta();
			for (auto& e: estimators) e.bootstrapTheta();
		_thetaOut.write("Bootstrap_" + toString(s + 1), "-", "-");

		logfile().endIndent();
	}

	// finish
	logfile().list("Total computation time for theta bootstrapping was ", timer.minutes());
	logfile().endIndent();
};

void TEstimateTheta::run() {
	_traverseBAMWindows();
	if (_genomeWide) {
		if (!_onlyBootstraps) {
			logfile().startIndent("Estimate theta based on a total of " + toString(_thetaEstimator.sizeWithData()) +
								  " sites:");
			_thetaEstimator.estimateTheta();
			for (auto& e: estimators) e.estimateTheta();
			// write estimates
			//std::string filename = _outputName + "_thetaGenomeWide.txt.gz";
			//_thetaOut.open(&_thetaEstimator, filename);
			if (_considerRegions) {
				_thetaOut.write("regions", "-", "-");
			} else {
				_thetaOut.write("genome-wide", "-", "-");
			}
			logfile().endIndent();
		}

		// bootstrap
		if (_numBootstraps > 0) { _bootstrapThetaEstimation(); }
	}
};

//-----------------------------------
// TEstimateThetaRatio
//-----------------------------------
TEstimateThetaRatio::TEstimateThetaRatio() : TGenome_windows(), _thetaEstimatorRatio() {
	// read the two regions to be used
	logfile().startIndent("Reading regions:");
	_initializeRegion(_region1, 1);
	_initializeRegion(_region2, 2);
};

void TEstimateThetaRatio::_initializeRegion(genometools::TBed &region, const int num) {
	logfile().startIndent((std::string) "Region " + std::to_string(num) + ":");
	std::string regionsFile = parameters().getParameter<std::string>("region" + std::to_string(num));
	logfile().listFlush("Reading regions ", num, " from file '", regionsFile, " (parameter 'region", num, "') ...");
	region.add(regionsFile, _bamFile.chromosomes());
	logfile().done();
	logfile().conclude("Read ", region.size(),  " sites on ", region.numChromosomesWithWindows(), " chromosomes.");
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
		UERROR("Failed to allocate sufficient memory to store the data for so many sites. Consider selecting fewer "
			   "regions or limiting to sites with a minimal depth (>=2 recommended).");
	}
	logfile().doneTime();
};

void TEstimateThetaRatio::run() {
	_traverseBAMWindows();
	_thetaEstimatorRatio.estimateRatio(_outputName);
};

}; // namespace GenomeTasks
