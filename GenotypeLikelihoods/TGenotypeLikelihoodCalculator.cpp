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
};

TGenotypeLikelihoodCalculator::TGenotypeLikelihoodCalculator(TParameters & params, const BAM::TReadGroups* ReadGroups, TLog* Logfile){
	_initialized = false;
	init(params, ReadGroups, Logfile);
};

void TGenotypeLikelihoodCalculator::init(TParameters & params, const BAM::TReadGroups* ReadGroups, TLog* Logfile){
	if(_initialized){
		throw "TGenotypeLikelihoodCalculator has already been initialized!";
	}

	//initialize PMD
	//--------------
	if(params.parameterExists("pmd")){
		std::vector<uint16_t> readGroupsWithoutDef;
		_pmdModels.initialize(params.getParameter<std::string>("pmd"), *ReadGroups, Logfile, readGroupsWithoutDef);

		//Warn if some read groups have no PMD definition
		if(readGroupsWithoutDef.size() > 0){
			Logfile->warning("The following read groups do not have PMD definitions: "
					         + concatenateString(ReadGroups->getNames(readGroupsWithoutDef), ", ")
							 + "!");
			if(!params.parameterExists("allowReadGroupsWithoutPMD")){
				throw "PMD is only defined for a subset of read groups. Did you use the wrong PMD file? (use allowReadGroupsWithoutPMD to ignore)";
			}
		}
	} else {
		Logfile->list("Assuming there is no PMD in the data. (use 'pmd' to add PMD definitions)");
	}

	//initialize sequencing errors
	//----------------------------
	if(params.parameterExists("recal")){
		_sequencingErrorModels.initializeFromFile(params.getParameter<std::string>("recal"), *ReadGroups, Logfile);

		//warn if some read groups have no recal definition
		std::vector<uint16_t> readGroupsWithoutRecal;
		std::vector<uint16_t> readGroupsLikelySingelEnd;
		_sequencingErrorModels.checkReadGroups(*ReadGroups, readGroupsWithoutRecal, readGroupsLikelySingelEnd);

		if(readGroupsWithoutRecal.size() > 0){
			Logfile->warning("The following read groups do not have recal definitions: "
					         + concatenateString(ReadGroups->getNames(readGroupsWithoutRecal), ", ")
					         + "!");
			if(!params.parameterExists("allowReadGroupsWithoutRecal")){
				throw "PMD is only defined for a subset of read groups. Did you use the wrong PMD file? (use allowReadGroupsWithoutRecal to ignore)";
			}
		}

		//Report if some read groups have only single-end definitions
		if(readGroupsLikelySingelEnd.size() > 0){
			Logfile->list("Read groups assumed single-end (no recal for second mate): "
					      + concatenateString(ReadGroups->getNames(readGroupsLikelySingelEnd), ", ")
						  + ".");
		}
	} else {
		Logfile->list("Assuming that error rates in BAM files are correct. (use 'recal' to add recalibration)");
	}

	_initialized = true;
};

bool TGenotypeLikelihoodCalculator::hasPMD() const{
	return _pmdModels.hasPMD();
};

bool TGenotypeLikelihoodCalculator::recalibrationChangesQualities() const{
	return _sequencingErrorModels.recalibrationChangesQualities();
};

BAM::ErrorRate TGenotypeLikelihoodCalculator::getErrorRate(const BAM::TSequencedBase & base) const{
	return _sequencingErrorModels.getErrorRate(base);
};

BAM::ErrorRate TGenotypeLikelihoodCalculator::getErrorWithPMD(const BAM::TSequencedBase & base) const{
	if(base.base == BAM::N){
		return 1.0;
	} else {
		//calculate base likelihoods with PMD
		static TBaseData tmpBaseData;
		_calculator.fillBaseLikelihoods(base, tmpBaseData, _pmdModels, _sequencingErrorModels);

		//return 1 - P(base|base) as in mapdamage2
		return 1.0 - tmpBaseData[base.base];
	}
};

BAM::PhredIntErrorRate TGenotypeLikelihoodCalculator::getPhredInt(const BAM::TSequencedBase & base) const{
	return _sequencingErrorModels.getPhredInt(base);
};

BAM::PhredIntErrorRate TGenotypeLikelihoodCalculator::getPhredIntWithPMD(const BAM::TSequencedBase & base) const{
	if(base.base == BAM::N){
		return BAM::PhredIntErrorRate::min();
	} else {
		return getErrorWithPMD(base);
	}
};

void TGenotypeLikelihoodCalculator::recalibrate(BAM::TSequencedBase & base) const{
	_sequencingErrorModels.recalibrate(base);
};

void TGenotypeLikelihoodCalculator::recalibrateWithPMD(BAM::TSequencedBase & base) const{
	base.recalibratedQualityAsPhredInt = getPhredIntWithPMD(base);
};

void TGenotypeLikelihoodCalculator::recalibrate(std::vector<BAM::TSequencedBase> & bases) const{
	_sequencingErrorModels.recalibrate(bases);
};

void TGenotypeLikelihoodCalculator::recalibrateWithPMD(std::vector<BAM::TSequencedBase> & bases) const{
	for(auto& b : bases){
		recalibrateWithPMD(b);
	}
};

double TGenotypeLikelihoodCalculator::calculateLogPMDS(const BAM::TSequencedBase & base, const BAM::Base & ref, const double & pi) const{
	//get base likelihoods
	static TBaseData baseLikelihoodsNoPMD;
	_sequencingErrorModels.fillBaseLikelihoods(base, baseLikelihoodsNoPMD);

	static TBaseData baseLikelihoods;
	_pmdModels.fillBaseLikelihoods(base, baseLikelihoodsNoPMD, baseLikelihoods);

	//calculate PMDS: true base in read == ref with prob. (1-pi) and different with prob. pi/3
	static TBaseData tmpBaseData;
	tmpBaseData.set(ref, pi);

	return log(baseLikelihoods.weightedSum(tmpBaseData) / baseLikelihoodsNoPMD.weightedSum(tmpBaseData));
};

void TGenotypeLikelihoodCalculator::calculateGenotypeLikelihoods(const TSite &site, TGenotypeLikelihoods &genotypeLikelihoods) const {
	_calculator.fillGenotypeLikelihoods(site, genotypeLikelihoods, _pmdModels, _sequencingErrorModels);
};

}; //end namespace

