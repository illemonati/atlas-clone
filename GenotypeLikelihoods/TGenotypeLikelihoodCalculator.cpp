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

TGenotypeLikelihoodCalculator::TGenotypeLikelihoodCalculator(const BAM::RGInfo::TReadGroupInfo& RGInfo) {
		_pmd.initialize(RGInfo);
		_recal.initialize(RGInfo);
};

coretools::Probability TGenotypeLikelihoodCalculator::errorWithPMD(const BAM::TSequencedBase &base) const {
	if (base.base == genometools::Base::N) { return coretools::Probability::highest(); }
	// calculate base likelihoods with PMD

	return _pmd.P_dij(base, _recal.P_dij(base))[base.base].complement();
};

genometools::PhredIntProbability TGenotypeLikelihoodCalculator::phredIntWithPMD(const BAM::TSequencedBase & base) const{
	if(base.base == genometools::Base::N){
		return genometools::PhredIntProbability::min();
	} else {
		return genometools::PhredIntProbability(errorWithPMD(base));
	}
};

void TGenotypeLikelihoodCalculator::recalibrateWithPMD(BAM::TSequencedBase & base) const{
	base.recalibratedQualityAsPhredInt = phredIntWithPMD(base);
};

void TGenotypeLikelihoodCalculator::recalibrateWithPMD(BAM::TAlignment& aln) const{
	for(auto& b : aln){
		recalibrateWithPMD(b);
	}
};

double TGenotypeLikelihoodCalculator::calculateLogPMDS(const BAM::TSequencedBase & base, const genometools::Base & ref, const coretools::Probability & pi) const{
	//get base likelihoods
	const TBaseLikelihoods baseLikelihoodsNoPMD = _recal.P_dij(base);

	const TBaseLikelihoods baseLikelihoods = _pmd.P_dij(base, baseLikelihoodsNoPMD);

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
			baseLikelihoods.push_back(_pmd.P_dij(d, _recal.P_dij(d)));
		}

		// calculate genotype likelihoods
		return getGLH(baseLikelihoods, site.depth());
};

}; //end namespace

