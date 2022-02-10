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

TGenotypeLikelihoodCalculator::TGenotypeLikelihoodCalculator(coretools::TParameters & params, const BAM::TReadGroups* ReadGroups, coretools::TLog* Logfile){
	_initialized = false;
	init(params, ReadGroups, Logfile);
};

void TGenotypeLikelihoodCalculator::init(coretools::TParameters & params, const BAM::TReadGroups* ReadGroups, coretools::TLog* Logfile){
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
					         + coretools::str::concatenateString(ReadGroups->getNames(readGroupsWithoutDef), ", ")
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
		std::string recalString = params.getParameter<std::string>("recal");

		//check if it is recal string
		if(_sequencingErrorModels.recalStringIsLikelyAModel(recalString)){
			//assume it is a model string
			Logfile->startIndent("Parsing common recal model for all read groups:");
			Logfile->list("Provided model (parameter 'recal'): " + recalString);

			std::string rhoString = params.getParameterWithDefault<std::string>("rho", "default");
			if(params.parameterExists("rho")){
				Logfile->list("Provided rho (parameter 'rho'): " + rhoString);
			} else {
				Logfile->list("Will use default rho. (use 'rho' to change)");
			}
			_sequencingErrorModels.initialize(recalString, rhoString, *ReadGroups, Logfile);
		} else {
			//assume it it a recal file
			Logfile->startIndent("Initializing recal models from file '" + recalString + "' (parameter 'recal'):");
			_sequencingErrorModels.initializeFromFile(recalString, *ReadGroups, Logfile);

			//warn if some read groups have no recal definition
			std::vector<uint16_t> readGroupsWithoutRecal;
			std::vector<uint16_t> readGroupsLikelySingelEnd;
			_sequencingErrorModels.checkReadGroups(*ReadGroups, readGroupsWithoutRecal, readGroupsLikelySingelEnd);

			if(readGroupsWithoutRecal.size() > 0){
				Logfile->warning("The following read groups do not have recal definitions: "
								 + coretools::str::concatenateString(ReadGroups->getNames(readGroupsWithoutRecal), ", ")
								 + "!");
				if(!params.parameterExists("allowReadGroupsWithoutRecal")){
					throw "Recal models are only defined for a subset of read groups. Did you use the wrong recal file? (use allowReadGroupsWithoutRecal to ignore)";
				}
			}

			//Report if some read groups have only single-end definitions
			if(readGroupsLikelySingelEnd.size() > 0){
				Logfile->list("Read groups assumed single-end (no recal for second mate): "
							  + coretools::str::concatenateString(ReadGroups->getNames(readGroupsLikelySingelEnd), ", ")
							  + ".");
			}
			Logfile->endIndent();
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

coretools::Probability TGenotypeLikelihoodCalculator::getErrorRate(const BAM::TSequencedBase & base) const{
	return _sequencingErrorModels.getErrorRate(base);
};

coretools::Probability TGenotypeLikelihoodCalculator::getErrorWithPMD(const BAM::TSequencedBase & base) const{
	if(base.base == genometools::Base::N){
		return Probability::highest();
	} else {
		//calculate base likelihoods with PMD
		static TBaseLikelihoods tmpBaseData;
		_calculator.fillBaseLikelihoods(base, tmpBaseData, _pmdModels, _sequencingErrorModels);

		//return 1 - P(base|base) as in mapdamage2
		return tmpBaseData[base.base].complement();
	}
};

genometools::PhredIntProbability TGenotypeLikelihoodCalculator::getPhredInt(const BAM::TSequencedBase & base) const{
	return _sequencingErrorModels.getPhredInt(base);
};

genometools::PhredIntProbability TGenotypeLikelihoodCalculator::getPhredIntWithPMD(const BAM::TSequencedBase & base) const{
	if(base.base == genometools::Base::N){
		return genometools::PhredIntProbability::min();
	} else {
		return genometools::PhredIntProbability(getErrorWithPMD(base));
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

double TGenotypeLikelihoodCalculator::calculateLogPMDS(const BAM::TSequencedBase & base, const genometools::Base & ref, const coretools::Probability & pi) const{
	//get base likelihoods
	static TBaseLikelihoods baseLikelihoodsNoPMD;
	_sequencingErrorModels.fillBaseLikelihoods(base, baseLikelihoodsNoPMD);

	static TBaseLikelihoods baseLikelihoods;
	_pmdModels.fillBaseLikelihoods(base, baseLikelihoodsNoPMD, baseLikelihoods);

	//calculate PMDS: true base in read == ref with prob. (1-pi) and different with prob. pi/3
	static TBaseLikelihoods tmpBaseData;
	tmpBaseData.setFromError(ref, pi);

	return log(baseLikelihoods.weightedSum(tmpBaseData) / baseLikelihoodsNoPMD.weightedSum(tmpBaseData));
};

void TGenotypeLikelihoodCalculator::calculateGenotypeLikelihoods(const TSite &site, TGenotypeLikelihoods &genotypeLikelihoods) const {
	_calculator.fillGenotypeLikelihoods(site, genotypeLikelihoods, _pmdModels, _sequencingErrorModels);
};

}; //end namespace

