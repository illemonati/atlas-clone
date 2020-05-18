/*
 * TGenotypeLikelihoods.cpp
 *
 *  Created on: May 8, 2020
 *      Author: phaentu
 */

#include <TGenotypeLikelihoodCalculator.h>

namespace GenotypeLikelihoods{

TGenotypeLikelihoodCalculator::TGenotypeLikelihoodCalculator(){
	initialized = false;
	logfile = nullptr;
	readGroups = nullptr;
	readGroupMap = nullptr;
};

TGenotypeLikelihoodCalculator::TGenotypeLikelihoodCalculator(TParameters & params, TReadGroups* ReadGroups, TLog* Logfile){
	initialized = false;
	init(params, ReadGroups, Logfile);
};

TGenotypeLikelihoodCalculator::~TGenotypeLikelihoodCalculator(){
	if(initialized){
		delete readGroupMap;
	}
};

void TGenotypeLikelihoodCalculator::init(TParameters & params, TReadGroups* ReadGroups, TLog* Logfile){
	if(initialized){
		throw "TGenotypeLikelihoodCalculator has alre<ady been initialized!";
	}
	logfile = Logfile;
	readGroups = ReadGroups;
	readGroupMap = new TReadGroupMap(ReadGroups);

	//initialize PMD
	pmd.initialize(params, *readGroups, logfile);

	//initialize sequencing errors
	if(params.parameterExists("recal")){
		sequencingErrorModels.createModels(params.getParameterString("recal"), readGroups, readGroupMap, logfile);
	} else {
		logfile->list("Assuming that error rates in BAM files are correct (no recalibration).");
	}

	//initialize storage to minimum size
	baseLikelihoods.resize(1);

	initialized = true;
};


double TGenotypeLikelihoodCalculator::getErrorRate(const TBaseData & base){
	return sequencingErrorModels.getErrorRate(base);
};

double TGenotypeLikelihoodCalculator::getErrorWithPMD(const TBaseData & base){
	if(base.base == N){
		return 1.0;
	} else {
		//calculate base likelihoods with PMD
		sequencingErrorModels.calculateBaseLikelihoods(base, baseLikelihoodsNoPMD);
		pmd.calculateBaseLikelihoods(base, baseLikelihoodsNoPMD, baseLikelihoods[0]);

		//return 1 - P(base|base) as in mapdamage2
		return 1.0 - baseLikelihoods[base.base];
	}
};

uint8_t TGenotypeLikelihoodCalculator::getPhredInt(const TBaseData & base){
	return sequencingErrorModels.getPhredInt(base);
};

uint8_t TGenotypeLikelihoodCalculator::getPhredIntWithPMD(const TBaseData & base){
	if(base.base == N){
		return 0;
	} else {
		return sequencingErrorModels.qualMap.errorToPhredInt(getErrorWithPMD(base));
	}
};

void TGenotypeLikelihoodCalculator::recalibrate(TBaseData & base){
	sequencingErrorModels.recalibrate(base);
};

void TGenotypeLikelihoodCalculator::recalibrateWithPMD(TBaseData & base){
	base.recalibratedQualityAsPhredInt = getPhredIntWithPMD(base);
};

void TGenotypeLikelihoodCalculator::recalibrate(TBase* bases, const uint16_t  length){
	sequencingErrorModels.recalibrate(bases, length);
};

void TGenotypeLikelihoodCalculator::recalibrateWithPMD(TBase* bases, const uint16_t  length){
	for(uint16_t i=0; i<length; ++i){
		recalibrateWithPMD(bases[i].data);
	}
};

void TGenotypeLikelihoodCalculator::calculateGenotypeLikelihoods(const std::vector<TBase*> bases, TGenotypeLikelihoods & genotypeLikelihoods){
	//ensure base likelihoods have proper size
	if(baseLikelihoods.size() < bases.size()){
		baseLikelihoods.resize(bases.size());
	}

	if(bases.empty()){
		genotypeLikelihoods.reset();
	} else {
		//calculate base likelihoods P(d|b, D, epsilon) = \sum_{\bar{b}} P(\bar{b}|b, D)P(d|\bar{b}, \epsilon)
		for(size_t i=0; i<bases.size(); ++i){
			sequencingErrorModels.calculateBaseLikelihoods(bases[i]->data, baseLikelihoodsNoPMD);
			pmd.calculateBaseLikelihoods(bases[i]->data, baseLikelihoodsNoPMD, baseLikelihoods[i]);
		}

		//calculate genotype likelihoods
		genotypeLikelihoods.fill(baseLikelihoods, bases.size());
	}
};




}; //end namespace

