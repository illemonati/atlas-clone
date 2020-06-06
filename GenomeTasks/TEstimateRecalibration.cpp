/*
 * TestimateRecalibration.cpp
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */


#include <TEstimateRecalibration.h>

namespace GenomeTasks{

//-----------------------------------------------------------
// TEstimateRecalibration_base
//-----------------------------------------------------------
TEstimateRecalibration_base::TEstimateRecalibration_base(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_windows(Parameters, Logfile, RandomGenerator){
	if(_genotypeLikelihoodCalculator.recalibrationChangesQualities() && !Parameters.parameterExists("rerecalibrate"))
		throw "Can not estimate recalibration: quality scores are already recalibrated while reading! (Use argument 'rerecalibrate' to overwrite this error)";

	//limit to sites with known alleles?
	if(Parameters.parameterExists("alleles")){
		_logfile->startIndent("Will limit analysis to sites with known alleles (parameter 'alleles'):");
		_openSiteSubset("alleles");
		_logfile->endIndent();
	} else {
		_logfile->list("Will use all sites without prior knowledge on alleles. (use 'alleles' to provide known alleles)");
	}

	//initialize maps
	if(Parameters.parameterExists("poolReadGroups")){
		std::string poolReadGroupsFile = Parameters.getParameterString("poolReadGroups");
		_logfile->startIndent("Will pool read groups (parameter 'poolReadGroups'):");
		_readGroupMap = std::make_unique<BAM::TReadGroupMap>(&(_bamFile.readGroups), poolReadGroupsFile, _logfile);
		_logfile->endIndent();
	} else  {
		_logfile->list("Will estimate recalibration parameters for each read group. (use 'poolReadGroups' to pool)");
		_readGroupMap = std::make_unique<BAM::TReadGroupMap>(&(_bamFile.readGroups));
	}

	//initialize recal estimator
	recalObjectEM = std::make_unique<GenotypeLikelihoods::TRecalibrationEMEstimator>(Parameters, _bamFile.readGroups, _logfile, _readGroupMap);
};

void TEstimateRecalibration_base::_handleWindow(){
	//add sites to recal estimator
	_window.estimateBaseFrequencies(_baseFreq);
	recalObjectEM->addNewWindow(_baseFreq);
	if(_subset){
		std::set<TSiteSubsetSite> thesePositions = _subset->getPositionInWindow(_window.startPos);
		for(auto& it : thesePositions){
			uint32_t internalPos = it.position - _window.startPos;
			if(_window[internalPos].hasData && it.ref==it.alt){
				recalObjectEM->addSite(_window[internalPos], _qualMap);
			}
		}
	} else {
		for(auto& s : _window){
			if(s.hasData){
				recalObjectEM->addSite(s, _qualMap);
			}
		}
	}
};


//-----------------------------------------------------------
// TEstimateRecalibration
//-----------------------------------------------------------
TEstimateRecalibration::TEstimateRecalibration(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TEstimateRecalibration_base(Parameters, Logfile, RandomGenerator){
	//write tmp tables?
	if(Parameters.parameterExists("writeTmpTables")){
		_writeTmpTables = true;
		_logfile->list("Will write intermediate estimates of EM and Newton-Raphson to file. (parameter 'writeTmpTables')");
	} else {
		_writeTmpTables = false;
		_logfile->list("Will not write intermediate estimates. (use 'writeTmpTables' to get this debug information)");
	}
};

void TEstimateRecalibration::estimateRecalibration(){
	//read data
	_traverseBAMWindows();

	//estimate recal parameters
	recalObjectEM->performEstimation(_outputName, _writeTmpTables);
};

//-----------------------------------------------------------
// TEstimateRecalibrationLL
//-----------------------------------------------------------
TEstimateRecalibrationLL::TEstimateRecalibrationLL(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TEstimateRecalibration_base(Parameters, Logfile, RandomGenerator){};

void TEstimateRecalibrationLL::calculateRecalibrationLL(){
	//read data
	_traverseBAMWindows();

	//estimate LL
	_logfile->list("LL = " + toString(recalObjectEM->calcLL()));
};

}; // end namespace
