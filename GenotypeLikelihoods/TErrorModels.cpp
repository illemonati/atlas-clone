/*
 * TGenotypeLikelihoods.cpp
 *
 *  Created on: May 8, 2020
 *      Author: phaentu
 */

#include "TErrorModels.h"
#include "GenotypeFunctions.h"
#include "TSite.h"
#include "coretools/algorithms.h"

namespace GenotypeLikelihoods{

TErrorModels::TErrorModels(const BAM::RGInfo::TReadGroupInfo& RGInfo) {
		_pmd.initialize(RGInfo);
		_recal.initialize(RGInfo);
};

coretools::Probability TErrorModels::errorWithPMD(const BAM::TSequencedBase &base) const {
	if (base.base == genometools::Base::N) { return coretools::Probability::highest(); }
	// calculate base likelihoods with PMD

	return _pmd.P_dij(base, _recal.P_dij(base))[base.base].complement();
};

genometools::PhredIntProbability TErrorModels::phredIntWithPMD(const BAM::TSequencedBase & base) const{
	if(base.base == genometools::Base::N){
		return genometools::PhredIntProbability::min();
	} else {
		return genometools::PhredIntProbability(errorWithPMD(base));
	}
};

void TErrorModels::recalibrateWithPMD(BAM::TSequencedBase & base) const{
	base.recalQuality = phredIntWithPMD(base);
};

void TErrorModels::recalibrateWithPMD(BAM::TAlignment& aln) const{
	for(auto& b : aln){
		recalibrateWithPMD(b);
	}
};

double TErrorModels::calculateLogPMDS(const BAM::TSequencedBase & base, const genometools::Base & ref, const coretools::Probability & pi) const{
	//get base likelihoods
	const TBaseLikelihoods baseLikelihoodsNoPMD = _recal.P_dij(base);

	const TBaseLikelihoods baseLikelihoods = _pmd.P_dij(base, baseLikelihoodsNoPMD);

	//calculate PMDS: true base in read == ref with prob. (1-pi) and different with prob. pi/3
	const TBaseLikelihoods tmpBaseData = fromError(ref, pi);

	return log(weightedSum(baseLikelihoods, tmpBaseData) / weightedSum(baseLikelihoodsNoPMD, tmpBaseData));
};

TGenotypeLikelihoods TErrorModels::calculateGenotypeLikelihoods(const TSite &site) const {
	if (site.empty()) { return TGenotypeLikelihoods{coretools::P(1.)}; }
	std::vector<TBaseLikelihoods> baseLikelihoods;
	baseLikelihoods.reserve(site.depth());
	// calculate base likelihoods P(d|b, D, epsilon) = \sum_{\bar{b}} P(\bar{b}|b, D)P(d|\bar{b}, \epsilon)
	for (const auto &d : site) { baseLikelihoods.push_back(_pmd.P_dij(d, _recal.P_dij(d))); }

	// calculate genotype likelihoods
	return getGLH(baseLikelihoods, site.depth());
};

}; //end namespace

