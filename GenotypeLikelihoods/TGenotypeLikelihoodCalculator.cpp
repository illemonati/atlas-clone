/*
 * TGenotypeLikelihoods.cpp
 *
 *  Created on: May 8, 2020
 *      Author: phaentu
 */

#include <TGenotypeLikelihoodCalculator.h>

namespace GenotypeLikelihoods{

TGenotypeLikelihoodCalculator::TGenotypeLikelihoodCalculator(){
	_initialized = false;
	_logfile = nullptr;
	_readGroups = nullptr;
	_readGroupMap = nullptr;
};

TGenotypeLikelihoodCalculator::TGenotypeLikelihoodCalculator(TParameters & params, BAM::TReadGroups* ReadGroups, TLog* Logfile){
	_initialized = false;
	init(params, ReadGroups, Logfile);
};

TGenotypeLikelihoodCalculator::~TGenotypeLikelihoodCalculator(){
	if(_initialized){
		delete _readGroupMap;
	}
};

void TGenotypeLikelihoodCalculator::init(TParameters & params, BAM::TReadGroups* ReadGroups, TLog* Logfile){
	if(_initialized){
		throw "TGenotypeLikelihoodCalculator has alre<ady been initialized!";
	}
	_logfile = Logfile;
	_readGroups = ReadGroups;
	_readGroupMap = new BAM::TReadGroupMap(ReadGroups);

	//initialize PMD
	_pmd.initialize(params, *_readGroups, _logfile);

	//initialize sequencing errors
	if(params.parameterExists("recal")){
		_sequencingErrorModels.createModels(params.getParameterString("recal"), _readGroups, _readGroupMap, _logfile);
	} else {
		_logfile->list("Assuming that error rates in BAM files are correct (no recalibration).");
	}

	//initialize storage to minimum size
	_baseLikelihoods.resize(1);

	_initialized = true;
};

bool TGenotypeLikelihoodCalculator::hasPMD() const{
	return _pmd.hasPMD();
};


bool TGenotypeLikelihoodCalculator::recalibrationChangesQualities() const{
	_sequencingErrorModels.recalibrationChangesQualities();
};

double TGenotypeLikelihoodCalculator::getErrorRate(const TBase & base) const{
	return _sequencingErrorModels.getErrorRate(base);
};

double TGenotypeLikelihoodCalculator::getErrorWithPMD(const TBase & base) const{
	if(base.base == N){
		return 1.0;
	} else {
		//calculate base likelihoods with PMD
		_sequencingErrorModels.calculateBaseLikelihoods(base, _baseLikelihoodsNoPMD);
		_pmd.calculateBaseLikelihoods(base, _baseLikelihoodsNoPMD, _baseLikelihoods[0]);

		//return 1 - P(base|base) as in mapdamage2
		return 1.0 - _baseLikelihoods[base.base];
	}
};

uint8_t TGenotypeLikelihoodCalculator::getPhredInt(const TBase & base) const{
	return _sequencingErrorModels.getPhredInt(base);
};

uint8_t TGenotypeLikelihoodCalculator::getPhredIntWithPMD(const TBase & base) const{
	if(base.base == N){
		return 0;
	} else {
		return _sequencingErrorModels.qualMap.errorToPhredInt(getErrorWithPMD(base));
	}
};

void TGenotypeLikelihoodCalculator::recalibrate(TBase & base) const{
	_sequencingErrorModels.recalibrate(base);
};

void TGenotypeLikelihoodCalculator::recalibrateWithPMD(TBase & base) const{
	base.recalibratedQualityAsPhredInt = getPhredIntWithPMD(base);
};

void TGenotypeLikelihoodCalculator::recalibrate(std::vector<TBase> & bases) const{
	_sequencingErrorModels.recalibrate(bases);
};

void TGenotypeLikelihoodCalculator::recalibrateWithPMD(std::vector<TBase> & bases) const{
	for(auto& b : bases){
		recalibrateWithPMD(b);
	}
};

void TGenotypeLikelihoodCalculator::calculateGenotypeLikelihoods(const std::vector<TBase*> bases, TGenotypeLikelihoods & genotypeLikelihoods) const{
	//ensure base likelihoods have proper size
	if(_baseLikelihoods.size() < bases.size()){
		_baseLikelihoods.resize(bases.size());
	}

	if(bases.empty()){
		genotypeLikelihoods.reset();
	} else {
		//calculate base likelihoods P(d|b, D, epsilon) = \sum_{\bar{b}} P(\bar{b}|b, D)P(d|\bar{b}, \epsilon)
		for(size_t i=0; i<bases.size(); ++i){
			_sequencingErrorModels.calculateBaseLikelihoods(*bases[i], _baseLikelihoodsNoPMD);
			_pmd.calculateBaseLikelihoods(*bases[i], _baseLikelihoodsNoPMD, _baseLikelihoods[i]);
		}

		//calculate genotype likelihoods
		genotypeLikelihoods.fill(_baseLikelihoods, bases.size());
	}
};

double TGenotypeLikelihoodCalculator::calculatePMDS(const TBase & base, const Base ref, const double pi) const{
	//get base likelihoods
	_sequencingErrorModels.calculateBaseLikelihoods(base, _baseLikelihoodsNoPMD);
	_pmd.calculateBaseLikelihoods(base, _baseLikelihoodsNoPMD, _baseLikelihoods[0]);

	//calculate PMDS: true base in read == ref with prob. (1-pi) and different with prob. pi/3
	_tmpBaseData.set(ref, pi);

	return log(_baseLikelihoods[0].weightedSum(_tmpBaseData) / _baseLikelihoodsNoPMD.weightedSum(_tmpBaseData));
};


}; //end namespace

