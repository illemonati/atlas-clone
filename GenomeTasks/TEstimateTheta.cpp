/*
 * TEstimateTheta.cpp
 *
 *  Created on: Jun 4, 2020
 *      Author: phaentu
 */

#include "TEstimateTheta.h"

namespace GenomeTasks{

//-----------------------------------
// TEstimateTheta_base
//-----------------------------------
TEstimateTheta_base::TEstimateTheta_base(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):
	TGenome_windows(Parameters, Logfile, RandomGenerator),
	_thetaEstimator(Parameters, Logfile, RandomGenerator){

};

void TEstimateTheta_base::_addSites(TWindow_base & window, GenotypeLikelihoods::TThetaEstimator & thetaEstimator){
	_logfile->listFlushTime("Calculating genotype likelihoods ...");
	for(auto& s : window){
		_genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(s._bases, _genoLik);
		thetaEstimator.add(s, _genoLik);
	}
	_logfile->doneTime();
};

void TEstimateTheta_base::_addSites(){
	_addSites(_window, _thetaEstimator);
};

//-----------------------------------
// TEstimateTheta
//-----------------------------------
TEstimateTheta::TEstimateTheta(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TEstimateTheta_base(Parameters, Logfile, RandomGenerator){
	//open output file
	std::string filename = _outputName + "_thetaPerWindow.txt.gz";
	_thetaOut.open(&_thetaEstimator, filename, _logfile);
	_logfile->list("Will write theta estimates to file '" + filename + "'.");

	//print all windows?
	if(Parameters.parameterExists("printAll")){
		_printAll = true;
		_logfile->list("Will print all windows, also those for which no estimation was possible. (parameter 'printAl')");
	} else {
		_printAll = false;
		_logfile->list("Will only print windows for which estimation was possible. (use 'printAll' to print all)");
	}
};

bool TEstimateTheta::_estimateTheta(){
	_addSites();

	//estimate Theta
	return _thetaEstimator.estimateTheta();
};

void TEstimateTheta::_handleWindow(){
	_logfile->startIndent("Estimating Theta:");
	if(_estimateTheta() || _printAll){
		_thetaOut.write(_window._chrName, _window.startPos, _window.endPos);
	}
	_thetaEstimator.clear();
	_logfile->endIndent();
};

void TEstimateTheta::estimateTheta(){
	_traverseBAMWindows();
};

//-----------------------------------
// TEstimateThetaGenomeWide
//-----------------------------------
TEstimateThetaGenomeWide::TEstimateThetaGenomeWide(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TEstimateTheta_base(Parameters, Logfile, RandomGenerator){
	if(_considerRegions){
		_logfile->list("Estimating theta at specific sites. (parameter 'regions')");
	} else {
		_logfile->list("Estimating theta genome-wide. (use 'regions' to limit)");
	}

	//bootstraps
	_numBootstraps = Parameters.getParameterIntWithDefault("bootstraps", 0);
	if(_numBootstraps > 0){
		_logfile->list("Will estimate theta fpr ", _numBootstraps, " bootstrap replicates. (parameter 'bootstraps')");
	} else {
		_logfile->list("Will not conduct any bootstrap replicates. (use 'bootstraps' to request)");
	}

	if(Parameters.parameterExists("onlyBootstrap")){
		_onlyBootstraps = true;
		_logfile->list("Will only ...");
	} else {
		_onlyBootstraps = true;
		_logfile->list("Will only ...");
	}
};

void TEstimateThetaGenomeWide::_handleWindow(){
	//add sites to data structure
	try{
		_addSites();
	} catch(...){
		throw "Failed to allocate sufficient memory to store the data for so many sites. Consider limiting the analysis (e.g. parameters 'regions'or 'chr') or limiting to sites with depth >=2 (parameter 'minDepth', recommended).";
	}
};

void TEstimateThetaGenomeWide::_bootstrapThetaEstimation(){
	_logfile->startIndent("Generating " + toString(_numBootstraps) + " bootstrap estimates of theta:");

	//measure runtime
	TTimer timer;

	//loop over bootstraps
	for(int s=0; s<_numBootstraps; ++s){
		_logfile->startIndent("Bootstrap " + toString(s+1) + " of " + toString(_numBootstraps) + ":");

		//run bootstrap
		_thetaEstimator.bootstrapTheta();
		_thetaOut.write("Bootstrap_" + toString(s+1), "-", "-");

		_logfile->endIndent();
	}

	//finish
	_logfile->list("Total computation time for theta bootstrapping was ", timer.minutes());
	_logfile->endIndent();
};

void TEstimateThetaGenomeWide::estimateThetaGenomeWide(){
	//read data
	_traverseBAMWindows();

	if(!_onlyBootstraps){
		_logfile->startIndent("Estimate theta based on a total of " + toString(_thetaEstimator.sizeWithData()) + " sites:");
		_thetaEstimator.estimateTheta();

		//write estimates
		std::string filename = _outputName + "_thetaGenomeWide.txt.gz";
		_thetaOut.open(&_thetaEstimator, filename, _logfile);
		_logfile->list("Will write theta estimates to file '" + filename + "'.");
		if(_considerRegions){
			_thetaOut.write("regions", "-", "-");
		} else {
			_thetaOut.write("genome-wide", "-", "-");
		}
		_logfile->endIndent();
	}

	//bootstrap
	if(_numBootstraps > 0){
		_bootstrapThetaEstimation();
	}
};

//-----------------------------------
// TEstimateThetaLLSurface
//-----------------------------------
TEstimateThetaLLSurface::TEstimateThetaLLSurface(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TEstimateTheta_base(Parameters, Logfile, RandomGenerator){
	_steps = Parameters.getParameterIntWithDefault("steps", 100);
	_logfile->list("Will calculate the LL-surface at ", _steps, " steps. (parameter 'steps')");
	if(_steps < 2){
		throw "Th enumber of steps must be >= 2!";
	}
};

void TEstimateThetaLLSurface::_handleWindow(){
	_logfile->startIndent("Calculating likelihood surface for Theta:");

	//adding sites to estimator
	_addSites();

	//open file
	std::string filename = _outputName + _window.chrName() + "_" + toString(_window.startPos) + "_LLsurface.txt";
	_logfile->listFlushTime("Writing LL surface to file '" + filename + "' ...");
	TOutputFile out(filename);

	_thetaEstimator.calcLikelihoodSurface(out, _steps);

	_thetaEstimator.clear();
	_logfile->doneTime();

	//finish
	_logfile->endIndent();
};

void TEstimateThetaLLSurface::estimateThetaLLSurface(){
	_traverseBAMWindows();
};

//-----------------------------------
// TEstimateThetaDownsamplingQC
//-----------------------------------
TEstimateThetaDownsamplingQC::TEstimateThetaDownsamplingQC(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TEstimateTheta_base(Parameters, Logfile, RandomGenerator){
	//parse downsampling parameters
	Parameters.fillParameterIntoProbabilityVectorWithDefault("prob", downSampleProbVector, ',', "1.0,0.5,0.2,0.1,0.05,0.02,0.01");

	//report probabilities
	_logfile->list("Will estimate theta after downsampling reads with probabilities " + concatenateString(downSampleProbVector, ", ") + ". (parameter 'prob')");

	//check if full data is to be used (i.e. if prob = 1.0 is specified)
	_printFullData = false;
	if(*downSampleProbVector.begin() == 1.0){
		_printFullData = true;
		downSampleProbVector.erase(downSampleProbVector.begin());
	}

	//create windows and estimators for downsamples
	if(downSampleProbVector.size() > 0){
		for(size_t i=0; i<downSampleProbVector.size(); ++i){
			estimators.emplace_back(_thetaEstimator);
		}
	}

	//open output
	std::string filename = _outputName + "_thetaQC.txt.gz";
	_logfile->list("Will write theta estimates to file '" + filename + "'.");
	if(_printFullData){
		_thetaOut.addEstimator(&_thetaEstimator, "p1.0_");
	}
	for(size_t i=0; i<downSampleProbVector.size(); ++i){
		//assemble prefix without lagging zeros
		std::string prefix = toString(downSampleProbVector[i]);
		int pos = prefix.size() - 1;
		while(prefix[pos] == '0') pos--;
		prefix.erase(pos + 1, prefix.size() - pos - 1);
		prefix = "p" + prefix + "_";

		//now add estimator ot output file
		_thetaOut.addEstimator(&estimators[i], prefix);
	}
	_thetaOut.open(filename, _logfile);
};

void TEstimateThetaDownsamplingQC::_handleWindow(){
	//estimate on full data
	if(_printFullData){
		_logfile->startIndent("Estimating Theta on full data:");
		_addSites();
		_thetaEstimator.estimateTheta();
		_logfile->endIndent();
	}

	for(size_t i=0; i<downSampleProbVector.size(); ++i){
		double p = downSampleProbVector[i];
		_logfile->startIndent("Estimating Theta on downsampled data (p = " + toString(p) + "):");

		//downsample
		_logfile->listFlush("Downsampling reads with p = " + toString(p) + " ...");
		if(_subset)
			destination.downsampleFromOther(_window, *_subset, _readUpToDepth, p, _randomGenerator);
		else
			destination.downsampleFromOther(_window, _readUpToDepth, p, _randomGenerator);
		_logfile->done();

		//apply filters
		_applyWindowFilters(destination);

		//and estimate
		_addSites(destination, estimators[i]);
		estimators[i].estimateTheta();
		_logfile->endIndent();
	}

	//write output
	_thetaOut.write(_window._chrName, _window.startPos, _window.endPos);

	//clear
	_thetaEstimator.clear();
	for(auto& e : estimators){
		e.clear();
	}
};

void TEstimateThetaDownsamplingQC::runQC(){
	_traverseBAMWindows();
};

//-----------------------------------
// TEstimateThetaRatio
//-----------------------------------
TEstimateThetaRatio::TEstimateThetaRatio(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_windows(Parameters, Logfile, RandomGenerator){
	//read the two regions to be used
	_logfile->startIndent("Reading regions:");
	_initializeRegion(Parameters, _region1, '1');
	_initializeRegion(Parameters, _region1, '2');
};

void TEstimateThetaRatio::_initializeRegion(TParameters & Parameters, BAM::TBedReaderWindows* region, const char num){
	_logfile->startIndent("Region " + num + ":");
	std::string regionsFile = Parameters.getParameterString("regions" + num);
	_logfile->list("Will read regions from file '" + regionsFile  + ". (parameter 'region" + num + "')");
	uint32_t siteLimit = 0;
	if(Parameters.parameterExists("siteLimit" + num)){
		siteLimit = Parameters.getParameterInt("siteLimit" + num);
		_logfile->list("Will only consider the first ", siteLimit, " sites. (parameter 'siteLimit" + num + "')");
		if(siteLimit < 10)
			throw "Site limit cannot be smaller than 10!";
	} else {
		_logfile->list("Will consider all sites in regions file. (use 'siteLimit" + num + "' to limit)");
	}
	region = new BAM::TBedReaderWindows(regionsFile, _windowSize, _bamFile._chromosomes, siteLimit,_logfile);
	_logfile->done();
};

void TEstimateThetaRatio::_addSites(GenotypeLikelihoods::TThetaEstimatorData & data, BAM::TBedReaderWindows & regions){
	if(regions.hasPositionsInWindow(_window.startPos)){
		GenotypeLikelihoods::TGenotypeLikelihoods genoLik;
		std::vector<uint32_t> thesePos = regions.getPositionInWindow(_window.startPos);
		for(auto& p : thesePos){
		uint32_t internalPos = p - _window.startPos;
			if(internalPos < _window.length){
				_genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(_window._sites[internalPos]._bases, genoLik);
				data.add(_window._sites[internalPos], genoLik);
			}
		}
	}
};

void TEstimateThetaRatio::_handleWindow(){
	_region1->setChr(_window._chrName);
	_region2->setChr(_window._chrName);

	//adding sites to estimator
	_logfile->listFlushTime("Calculating genotype likelihoods ...");
	try{
		_addSites(*_thetaEstimatorRatio.pointerToDataContainer(), *_region1);
		_addSites(*_thetaEstimatorRatio.pointerToDataContainer2(), *_region2);
	} catch(...){
		throw "Failed to allocate sufficient memory to store the data for so many sites. Consider selecting fewer regions or limiting to sites with a minimal depth (>=2 recommended).";
	}
	_logfile->doneTime();
};

void TEstimateThetaRatio::estimateThetaRation(){
	_traverseBAMWindows();
	_thetaEstimatorRatio.estimateRatio(_outputName);
};

}; // end namespace
