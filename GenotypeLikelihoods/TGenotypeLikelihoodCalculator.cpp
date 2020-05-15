/*
 * TGenotypeLikelihoods.cpp
 *
 *  Created on: May 8, 2020
 *      Author: phaentu
 */

#include <TGenotypeLikelihoodCalculator.h>

namespace GenotypeLikelihoods{

TGenotypeLikelihoodCalculator::TGenotypeLikelihoodCalculator(TParameters & params, TReadGroups* ReadGroups, TLog* Logfile){
	logfile = Logfile;
	readGroups = ReadGroups;
	readGroupMap = new TReadGroupMap(ReadGroups);

	//initialize PMD
	pmd.initialize(params, *readGroups, logfile);

	//initialize sequencing errors
	if(params.parameterExists("recal")){
		errorModels.createModels(params.getParameterString("recal"), readGroups, readGroupMap, logfile);
	} else {
		logfile->list("Assuming that error rates in BAM files are correct (no recalibration).");
	}
};


void TGenotypeLikelihoodCalculator::calculateGenotypeLikelihoods(const std::vector<TBase*> bases, TGenotypeLikelihoods & genotypeLikelihoods){
	//ensure base likelihoods have proper size
	if(baseLikelihoods.size() < bases.size()){
		baseLikelihoods.resize(bases.size());
	}

	//calculate base likelihoods P(d|b, D, epsilon) = \sum_{\bar{b}} P(\bar{b}|b, D)P(d|\bar{b}, \epsilon)
	for(size_t i=0; i<bases.size(); ++i){
		errorModels.calculateBaseLikelihoods(bases[i]->data, baseLikelihoodsNoPMD);
		pmd.calculateBaseLikelihoods(bases[i]->data, baseLikelihoodsNoPMD, baseLikelihoods[i]);
	}

	//calculate genotype likelihoods
	genotypeLikelihoods.fill(baseLikelihoods, bases.size());
};




}; //end namespace

