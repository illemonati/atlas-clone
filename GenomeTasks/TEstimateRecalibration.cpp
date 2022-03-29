/*
 * TestimateRecalibration.cpp
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */


#include "TEstimateRecalibration.h"

namespace GenomeTasks{

//-----------------------------------------------------------
// TEstimateRecalibration_base
//-----------------------------------------------------------
TEstimateRecalibration::TEstimateRecalibration(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator):TGenome_windows(Parameters, Logfile, RandomGenerator){
	if(_genotypeLikelihoodCalculator.recalibrationChangesQualities() && !Parameters.parameterExists("rerecalibrate"))
		throw "Can not estimate recalibration: quality scores are already recalibrated while reading! (Use argument 'rerecalibrate' to overwrite this error)";

	//limit to sites with known alleles?
	if(Parameters.parameterExists("alleles")){
		_logfile->startIndent("Will limit analysis to sites with known genotypes (parameter 'genotypes'):");
		_openSiteSubset("alleles");
		_logfile->endIndent();
	} else {
		_logfile->list("Will use all sites without prior knowledge on alleles. (use 'genotypes' to provide known genotypes)");
	}

	//initialize maps
	if(Parameters.parameterExists("poolReadGroups")){
		std::string poolReadGroupsFile = Parameters.getParameter<std::string>("poolReadGroups");
		_logfile->startIndent("Will pool read groups (parameter 'poolReadGroups'):");
		_readGroupMap = new BAM::TReadGroupMap(_bamFile.readGroups(), poolReadGroupsFile, _logfile);
		_logfile->endIndent();
	} else  {
		_logfile->list("Will estimate recalibration parameters for each read group. (use 'poolReadGroups' to pool)");
		_readGroupMap = new BAM::TReadGroupMap(_bamFile.readGroups());
	}

	//initialize recal estimator
	recalObjectEM = std::make_unique<GenotypeLikelihoods::SequencingError::TRecalibrationEMEstimator>(&(_bamFile.readGroups()), _readGroupMap);
};

TEstimateRecalibration::~TEstimateRecalibration(){
	delete _readGroupMap;
};

void TEstimateRecalibration::_handleWindow(){
	//add sites to recal estimator
	if(_subset){
		std::set<GenotypeLikelihoods::TSiteSubsetSite> thesePositions = _subset->getPositionInWindow(_window);
		for(auto& it : thesePositions){
			uint32_t internalPos = it - _window.from();
			if(!_window[internalPos].empty() && it.ref() == it.alt()){
				recalObjectEM->addSite(_window[internalPos]);
			}
		}
	} else {
		for(auto& s : _window){
			if(!s.empty()){
				recalObjectEM->addSite(s);
			}
		}
	}
};

void TEstimateRecalibration::estimateRecalibration(){
	//read data
	_traverseBAMWindows();

	//estimate recal parameters
	recalObjectEM->performEstimation(_outputName, _genotypeLikelihoodCalculator.getSequencingErrorModelsMutable(), _genotypeLikelihoodCalculator.getPostMortemDamageModelsMutable());
};


}; // end namespace
