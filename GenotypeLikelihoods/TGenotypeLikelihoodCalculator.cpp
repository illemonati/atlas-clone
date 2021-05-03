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
		throw "TGenotypeLikelihoodCalculator has already been initialized!";
	}
	_logfile = Logfile;
	_readGroups = ReadGroups;

	//initialize storage to minimum size
	_baseLikelihoods.resize(1);

	//initialize PMD
	//--------------
	if(params.parameterExists("pmd")){
		std::vector<uint16_t> readGroupsWithoutDef;
		_pmd.initialize(params.getParameterString("pmd"), *_readGroups, _logfile, readGroupsWithoutDef);

		//Warn if some read groups have no PMD definition
		if(readGroupsWithoutDef.size() > 0){
			_logfile->warning("The following read groups do not have PMD definitions: "
					         + concatenateString(_readGroups->getNames(readGroupsWithoutDef), ", ")
							 + "!");
			if(!params.parameterExists("allowReadGroupsWithoutPMD")){
				throw "PMD is only defined for a subset of read groups. Did you use the wrong PMD file? (use allowReadGroupsWithoutPMD to ignore)";
			}
		}
	} else {
		_logfile->list("Assuming there is no PMD in the data. (use 'pmd' to add PMD definitions)");
	}

	//initialize sequencing errors
	//----------------------------
	if(params.parameterExists("recal")){
		std::vector<uint16_t> readGroupsWithoutRecal;
		std::vector<uint16_t> readGroupsLikelySingelEnd;
		_sequencingErrorModels.initializeFromFile(params.getParameterString("recal"), *_readGroups, _logfile, readGroupsWithoutRecal, readGroupsLikelySingelEnd);

		//warn if some read groups have no PMD definition
		if(readGroupsWithoutRecal.size() > 0){
			_logfile->warning("The following read groups do not have recal definitions: "
					         + concatenateString(_readGroups->getNames(readGroupsWithoutRecal), ", ")
					         + "!");
			if(!params.parameterExists("allowReadGroupsWithoutRecal")){
				throw "PMD is only defined for a subset of read groups. Did you use the wrong PMD file? (use allowReadGroupsWithoutRecal to ignore)";
			}
		}

		//Report if some read groups have only single-end definitions
		if(readGroupsLikelySingelEnd.size() > 0){
			_logfile->list("Read groups assumed single-end (no recal for second mate): "
					      + concatenateString(_readGroups->getNames(readGroupsLikelySingelEnd), ", ")
						  + ".");
		}
	} else if(){

	} else {
		_logfile->list("Assuming that error rates in BAM files are correct. (use 'recal' to add recalibration)");
	}

	_initialized = true;
};

bool TGenotypeLikelihoodCalculator::hasPMD() const{
	return _pmd.hasPMD();
};

bool TGenotypeLikelihoodCalculator::recalibrationChangesQualities() const{
	return _sequencingErrorModels.recalibrationChangesQualities();
};

double TGenotypeLikelihoodCalculator::getErrorRate(const BAM::TBase & base) const{
	return _sequencingErrorModels.getErrorRate(base);
};

double TGenotypeLikelihoodCalculator::getErrorWithPMD(const BAM::TBase & base) const{
	if(base.base == N){
		return 1.0;
	} else {
		//calculate base likelihoods with PMD
		_sequencingErrorModels.calculateBaseLikelihoods(base, _baseLikelihoodsNoPMD);
		_pmd.calculateBaseLikelihoods(base, _baseLikelihoodsNoPMD, _baseLikelihoods[0]);

		//return 1 - P(base|base) as in mapdamage2
		return 1.0 - _baseLikelihoods[base.base][base.base];
	}
};

uint8_t TGenotypeLikelihoodCalculator::getPhredInt(const BAM::TBase & base) const{
	return _sequencingErrorModels.getPhredInt(base);
};

uint8_t TGenotypeLikelihoodCalculator::getPhredIntWithPMD(const BAM::TBase & base) const{
	if(base.base == N){
		return 0;
	} else {
		return _sequencingErrorModels.qualMap.errorToPhredInt(getErrorWithPMD(base));
	}
};

void TGenotypeLikelihoodCalculator::recalibrate(BAM::TBase & base) const{
	_sequencingErrorModels.recalibrate(base);
};

void TGenotypeLikelihoodCalculator::recalibrateWithPMD(BAM::TBase & base) const{
	base.recalibratedQualityAsPhredInt = getPhredIntWithPMD(base);
};

void TGenotypeLikelihoodCalculator::recalibrate(std::vector<BAM::TBase> & bases) const{
	_sequencingErrorModels.recalibrate(bases);
};

void TGenotypeLikelihoodCalculator::recalibrateWithPMD(std::vector<BAM::TBase> & bases) const{
	for(auto& b : bases){
		recalibrateWithPMD(b);
	}
};

double TGenotypeLikelihoodCalculator::calculateLogPMDS(const BAM::TBase & base, const Base ref, const double pi) const{
	//get base likelihoods
	_sequencingErrorModels.calculateBaseLikelihoods(base, _baseLikelihoodsNoPMD);
	_pmd.calculateBaseLikelihoods(base, _baseLikelihoodsNoPMD, _baseLikelihoods[0]);

	//calculate PMDS: true base in read == ref with prob. (1-pi) and different with prob. pi/3
	_tmpBaseData.set(ref, pi);

	return log(_baseLikelihoods[0].weightedSum(_tmpBaseData) / _baseLikelihoodsNoPMD.weightedSum(_tmpBaseData));
};

void TGenotypeLikelihoodCalculator::calculateGenotypeLikelihoods(const TSite &site, TGenotypeLikelihoods &genotypeLikelihoods) const {
    //ensure base likelihoods have proper size
    if(_baseLikelihoods.size() < site.depth()){
        _baseLikelihoods.resize(site.depth());
    }

    if(site.empty()){
        genotypeLikelihoods.reset();
    } else {
        //calculate base likelihoods P(d|b, D, epsilon) = \sum_{\bar{b}} P(\bar{b}|b, D)P(d|\bar{b}, \epsilon)
        for(size_t i=0; i<site.depth(); ++i){
            _sequencingErrorModels.calculateBaseLikelihoods(site[i], _baseLikelihoodsNoPMD);
            _pmd.calculateBaseLikelihoods(site[i], _baseLikelihoodsNoPMD, _baseLikelihoods[i]);
        }

        //calculate genotype likelihoods
        genotypeLikelihoods.fill(_baseLikelihoods, site.depth());
    }
}

/*
void TGenotypeLikelihoodCalculator::calculateGenotypeLikelihoods(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) const{
	//ensure base likelihoods have proper size
	if(_baseLikelihoods.size() < site.depth()){
		_baseLikelihoods.resize(site.depth());
	}

	if(site.empty()){
		genotypeLikelihoods.reset();
	} else {
		//calculate base likelihoods P(d|b, D, epsilon) = \sum_{\bar{b}} P(\bar{b}|b, D)P(d|\bar{b}, \epsilon)
		for(size_t i=0; i<site.depth(); ++i){
			_sequencingErrorModels.calculateBaseLikelihoods(site.at(i), _baseLikelihoodsNoPMD);
			_pmd.calculateBaseLikelihoods(site.at(i), _baseLikelihoodsNoPMD, _baseLikelihoods[i]);
		}

		//calculate genotype likelihoods
		genotypeLikelihoods.fill(_baseLikelihoods, site.depth());
	}
};
*/



}; //end namespace

