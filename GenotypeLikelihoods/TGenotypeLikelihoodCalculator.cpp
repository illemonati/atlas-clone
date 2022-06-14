/*
 * TGenotypeLikelihoods.cpp
 *
 *  Created on: May 8, 2020
 *      Author: phaentu
 */

#include "TGenotypeLikelihoodCalculator.h"

#include <math.h>
#include <stdint.h>
#include <exception>
#include <string>

#include "TGenotypeData.h"
#include "TLog.h"
#include "TParameters.h"
#include "TReadGroups.h"
#include "TSequencedBase.h"
#include "stringFunctions.h"

namespace GenotypeLikelihoods{

using coretools::instances::parameters;
using coretools::instances::logfile;

namespace impl {

bool isLikelyAModel(const std::string &RecalString) noexcept {
	// check if it contains a ';', ':', '[' or ']'
	return coretools::str::stringContainsAny(RecalString, ";:[]");
}
} // namespace impl

TGenotypeLikelihoodCalculator::TGenotypeLikelihoodCalculator(){
	_initialized = false;
};

TGenotypeLikelihoodCalculator::TGenotypeLikelihoodCalculator(const BAM::TReadGroups* ReadGroups){
	_initialized = false;
	init(ReadGroups);
};

void TGenotypeLikelihoodCalculator::init(const BAM::TReadGroups* ReadGroups){
	if(_initialized){
		throw "TGenotypeLikelihoodCalculator has already been initialized!";
	}

	//initialize PMD
	//--------------
	if(parameters().parameterExists("pmd")){
		std::vector<uint16_t> readGroupsWithoutDef;
		_pmdModels.initialize(parameters().getParameter<std::string>("pmd"), *ReadGroups, readGroupsWithoutDef);

		//Warn if some read groups have no PMD definition
		if(readGroupsWithoutDef.size() > 0){
			logfile().warning("The following read groups do not have PMD definitions: "
					         + coretools::str::concatenateString(ReadGroups->getNames(readGroupsWithoutDef), ", ")
							 + "!");
			if(!parameters().parameterExists("allowReadGroupsWithoutPMD")){
				throw "PMD is only defined for a subset of read groups. Did you use the wrong PMD file? (use allowReadGroupsWithoutPMD to ignore)";
			}
		}
	} else {
		logfile().list("Assuming there is no PMD in the data. (use 'pmd' to add PMD definitions)");
	}

	//initialize sequencing errors
	//----------------------------
	if(parameters().parameterExists("recal")){
		std::string recalString = parameters().getParameter<std::string>("recal");

		//check if it is recal string
		if(impl::isLikelyAModel(recalString)){
			//assume it is a model string
			logfile().startIndent("Parsing common recal model for all read groups:");
			logfile().list("Provided model (parameter 'recal'): " + recalString);

			std::string rhoString = parameters().getParameterWithDefault<std::string>("rho", "default");
			if(parameters().parameterExists("rho")){
				logfile().list("Provided rho (parameter 'rho'): " + rhoString);
			} else {
				logfile().list("Will use default rho. (use 'rho' to change)");
			}
			_sequencingErrorModels.initialize(recalString, rhoString, *ReadGroups);
		} else {
			//assume it it a recal file
			logfile().startIndent("Initializing recal models from file '" + recalString + "' (parameter 'recal'):");
			_sequencingErrorModels.initializeFromFile(recalString, *ReadGroups);
			//warn if some read groups have no recal definition
			std::vector<uint16_t> readGroupsWithoutRecal;
			std::vector<uint16_t> readGroupsLikelySingelEnd;

			_sequencingErrorModels.checkReadGroups(*ReadGroups, readGroupsWithoutRecal, readGroupsLikelySingelEnd);

			if(readGroupsWithoutRecal.size() > 0){
				logfile().warning("The following read groups do not have recal definitions: "
								 + coretools::str::concatenateString(ReadGroups->getNames(readGroupsWithoutRecal), ", ")
								 + "!");
				if(!parameters().parameterExists("allowReadGroupsWithoutRecal")){
					throw "Recal models are only defined for a subset of read groups. Did you use the wrong recal file? (use allowReadGroupsWithoutRecal to ignore)";
				}
			}

			//Report if some read groups have only single-end definitions
			if(readGroupsLikelySingelEnd.size() > 0){
				logfile().list("Read groups assumed single-end (no recal for second mate): "
							  + coretools::str::concatenateString(ReadGroups->getNames(readGroupsLikelySingelEnd), ", ")
							  + ".");
			}
			logfile().endIndent();
		}
	} else {
		logfile().list("Assuming that error rates in BAM files are correct. (use 'recal' to add recalibration)");
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

coretools::Probability TGenotypeLikelihoodCalculator::getErrorWithPMD(const BAM::TSequencedBase &base) const {
	if (base.base == genometools::Base::N) { return coretools::Probability::highest(); }
	// calculate base likelihoods with PMD
	TBaseLikelihoods tmpBaseData = _calculator.getBaseLikelihoods(base, _pmdModels, _sequencingErrorModels);

	// return 1 - P(base|base) as in mapdamage2
	return tmpBaseData[base.base].complement();
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
	const TBaseLikelihoods baseLikelihoodsNoPMD = _sequencingErrorModels.getBaseLikelihoods(base);

	const TBaseLikelihoods baseLikelihoods = _pmdModels.getBaseLikelihoods(base, baseLikelihoodsNoPMD);

	//calculate PMDS: true base in read == ref with prob. (1-pi) and different with prob. pi/3
	const TBaseLikelihoods tmpBaseData = fromError(ref, pi);

	return log(weightedSum(baseLikelihoods, tmpBaseData) / weightedSum(baseLikelihoodsNoPMD, tmpBaseData));
};

TGenotypeLikelihoods TGenotypeLikelihoodCalculator::calculateGenotypeLikelihoods(const TSite &site) const {
	return _calculator.getGenotypeLikelihoods(site, _pmdModels, _sequencingErrorModels);
};

}; //end namespace

