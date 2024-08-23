/*
 * TGenotypeLikelihoods.cpp
 *
 *  Created on: May 8, 2020
 *      Author: phaentu
 */

#include "TErrorModels.h"

#include "GenotypeFunctions.h"
#include "TSite.h"
#include "coretools/Main/TLog.h"
#include "coretools/algorithms.h"
#include "genometools/GenotypeContainers.h"
#include "genometools/GenotypeTypes.h"

namespace GenotypeLikelihoods{
using coretools::instances::logfile;
using genometools::TBaseLikelihoods;
using genometools::TGenotypeLikelihoods;

TErrorModels::TErrorModels(const BAM::RGInfo::TReadGroupInfo& RGInfo) {
		_pmd.initialize(RGInfo);
		_recal.initialize(RGInfo);
		logfile().startIndent("Using the following error Models:");
		for (size_t rg = 0; rg < RGInfo.size(); ++rg) {
			logfile().startIndent("Readgroup ", rg, ":");
			_recal.log(rg);
			_pmd.log(rg);
			logfile().endIndent();
		}
		logfile().endIndent();
};

coretools::Probability TErrorModels::errorWithPMD(const BAM::TSequencedBase &base) const {
	if (base.base == genometools::Base::N) { return coretools::Probability::highest(); }
	// calculate base likelihoods with PMD

	return _pmd.P_dij(base, _recal.P_dij(base))[base.base].complement();
};

coretools::PhredInt TErrorModels::phredIntWithPMD(const BAM::TSequencedBase & base) const{
	if(base.base == genometools::Base::N){
		return coretools::PhredInt::highest();
	} else {
		return coretools::PhredInt(errorWithPMD(base));
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
	using coretools::P;
	using genometools::Base;
	using genometools::Genotype;
		
	if (site.empty()) { return TGenotypeLikelihoods{coretools::P(1.)}; }
	std::vector<TBaseLikelihoods> baseLikelihoods;
	baseLikelihoods.reserve(site.depth());
	// calculate base likelihoods P(d|b, D, epsilon) = \sum_{\bar{b}} P(\bar{b}|b, D)P(d|\bar{b}, \epsilon)
	TGenotypeLikelihoods ret{P(1.)};

	// Normalize to max = 1
	double max = 1.;
	for (const auto &d : site) {
		const double max_inv = 1./max;
		max = 0.;
		const auto bLikes = _pmd.P_dij(d, _recal.P_dij(d));
		for (auto k = Base::min; k < Base::max; ++k) {
			const auto kk = genometools::genotype(k, k);
			ret[kk] *= P(bLikes[k]*max_inv);
			max = std::max<double>(ret[kk], max);
		}
		for (const auto kl: {Genotype::AC, Genotype::AG, Genotype::AT, Genotype::CG, Genotype::CT, Genotype::GT}) {
			const auto k = genometools::first(kl);
			const auto l = genometools::second(kl);
			ret[kl] *= P(0.5*(bLikes[k] + bLikes[l])*max_inv);
			max = std::max<double>(ret[kl], max);
		}
	}
	const double max_inv = 1./max;
	for (auto& r: ret) {
		r = P(r*max_inv);
	}

	return ret;

};

}; //end namespace
