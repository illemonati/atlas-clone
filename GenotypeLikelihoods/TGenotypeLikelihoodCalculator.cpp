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
#include "TReadGroupInfo.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "TReadGroups.h"
#include "TSequencedBase.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenotypeLikelihoods{

using coretools::instances::parameters;
using coretools::instances::logfile;

TGenotypeLikelihoodCalculator::TGenotypeLikelihoodCalculator(const BAM::TReadGroups* ReadGroups){
	//initialize PMD
	//--------------
	if(parameters().parameterExists("pmd")){
		std::vector<size_t> readGroupsWithoutDef = _pmdModels.initialize(parameters().getParameter<std::string>("pmd"), *ReadGroups);

		//Warn if some read groups have no PMD definition
		if(readGroupsWithoutDef.size() > 0){
			logfile().warning("The following read groups do not have PMD definitions: ",
					         coretools::str::concatenateString(ReadGroups->getNames(readGroupsWithoutDef), ", ")
							 + "!");
			if(!parameters().parameterExists("allowReadGroupsWithoutPMD")){
				UERROR("PMD is only defined for a subset of read groups. Did you use the wrong PMD file? (use allowReadGroupsWithoutPMD to ignore)");
			}
		}
	} else {
		logfile().list("Assuming there is no PMD in the data. (use 'pmd' to add PMD definitions)");
	}

	//initialize sequencing errors
	//----------------------------
	if (parameters().parameterExists(BAM::RGInfo::TReadGroupInfo::RGInfoArgument)) {
		BAM::RGInfo::TReadGroupInfo RGinfo(*ReadGroups, parameters().getParameter(BAM::RGInfo::TReadGroupInfo::RGInfoArgument));
		_sequencingErrorModels.initialize(RGinfo);
	} else if(parameters().parameterExists("recal")){
		std::string recalString = parameters().getParameter<std::string>("recal");

		//check if it is recal string

		if (std::filesystem::exists(recalString)) {
			UERROR("Initializing recal with recal-file is not supported anymore. Use Readgroup Info File in json-format instead!");
		}
		else {
			//assume it is a model string
			logfile().startIndent("Parsing common recal model for all read groups:");
			logfile().list("Provided model (parameter 'recal'): " + recalString);

			std::string rhoString = parameters().getParameterWithDefault<std::string>("rho", "default");
			if(parameters().parameterExists("rho")){
				logfile().list("Provided rho (parameter 'rho'): " + rhoString);
			} else {
				logfile().list("Will use default rho. (use 'rho' to change)");
			}
			_sequencingErrorModels.initialize(recalString, rhoString, ReadGroups->size());
		}
	} else {
		_sequencingErrorModels.initializeNoRecal(ReadGroups->size());
		logfile().list("Assuming that error rates in BAM files are correct. (use 'recal' to add recalibration)");
	}
};

bool TGenotypeLikelihoodCalculator::hasPMD() const{
	return _pmdModels.hasPMD();
};

bool TGenotypeLikelihoodCalculator::recalibrationChangesQualities() const{
	return _sequencingErrorModels.recalibrates();
};

coretools::Probability TGenotypeLikelihoodCalculator::errorRate(const BAM::TSequencedBase & base) const{
	return _sequencingErrorModels.errorRate(base);
};

coretools::Probability TGenotypeLikelihoodCalculator::errorWithPMD(const BAM::TSequencedBase &base) const {
	if (base.base == genometools::Base::N) { return coretools::Probability::highest(); }
	// calculate base likelihoods with PMD

	return _pmdModels.P_dij(base, _sequencingErrorModels.baseLikelihoods(base))[base.base].complement();
};

genometools::PhredIntProbability TGenotypeLikelihoodCalculator::phredInt(const BAM::TSequencedBase & base) const{
	return _sequencingErrorModels.phredInt(base);
};

genometools::PhredIntProbability TGenotypeLikelihoodCalculator::phredIntWithPMD(const BAM::TSequencedBase & base) const{
	if(base.base == genometools::Base::N){
		return genometools::PhredIntProbability::min();
	} else {
		return genometools::PhredIntProbability(errorWithPMD(base));
	}
};

void TGenotypeLikelihoodCalculator::recalibrate(BAM::TSequencedBase & base) const{
	_sequencingErrorModels.recalibrate(base);
};

void TGenotypeLikelihoodCalculator::recalibrateWithPMD(BAM::TSequencedBase & base) const{
	base.recalibratedQualityAsPhredInt = phredIntWithPMD(base);
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
	const TBaseLikelihoods baseLikelihoodsNoPMD = _sequencingErrorModels.baseLikelihoods(base);

	const TBaseLikelihoods baseLikelihoods = _pmdModels.P_dij(base, baseLikelihoodsNoPMD);

	//calculate PMDS: true base in read == ref with prob. (1-pi) and different with prob. pi/3
	const TBaseLikelihoods tmpBaseData = fromError(ref, pi);

	return log(weightedSum(baseLikelihoods, tmpBaseData) / weightedSum(baseLikelihoodsNoPMD, tmpBaseData));
};

TGenotypeLikelihoods TGenotypeLikelihoodCalculator::calculateGenotypeLikelihoods(const TSite &site) const {
		if (site.empty()) { return TGenotypeLikelihoods{1.}; }
		std::vector<TBaseLikelihoods> baseLikelihoods;
		baseLikelihoods.reserve(site.depth());
		// calculate base likelihoods P(d|b, D, epsilon) = \sum_{\bar{b}} P(\bar{b}|b, D)P(d|\bar{b}, \epsilon)
		for (const auto &d : site) {
			baseLikelihoods.push_back(_pmdModels.P_dij(d, _sequencingErrorModels.baseLikelihoods(d)));
		}

		// calculate genotype likelihoods
		return getGLH(baseLikelihoods, site.depth());
};

}; //end namespace

