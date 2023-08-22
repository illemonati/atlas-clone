/*
 * PMD/TModels.h
 *
 */

#ifndef GENOTYPELIKELIHOODS_PMD_MODEL_H_
#define GENOTYPELIKELIHOODS_PMD_MODEL_H_


#include "TGenotypeData.h"
namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }

namespace GenotypeLikelihoods::PMD {

struct TModel {
	virtual TBaseLikelihoods P_dij(const BAM::TSequencedBase &data,
								   const TBaseLikelihoods &P_dij_bbar) const noexcept    = 0;
	virtual TBaseProbabilities P_bbar(genometools::Base b, const BAM::TSequencedBase &data,
									  const TBaseLikelihoods &P_dij_bbar) const noexcept = 0;
	virtual TBaseProbabilities P_bbar(genometools::Genotype g, const BAM::TSequencedBase &data,
									  const TBaseLikelihoods &P_dij_bbar) const noexcept = 0;
};

struct TNoPMD final : public TModel {
	TBaseLikelihoods P_dij(const BAM::TSequencedBase &, const TBaseLikelihoods &P_dij_bbar) const noexcept override {
		return P_dij_bbar;
	}
	TBaseProbabilities P_bbar(genometools::Base b, const BAM::TSequencedBase &data,
							  const TBaseLikelihoods &P_dij_bbar) const noexcept override;
	TBaseProbabilities P_bbar(genometools::Genotype g, const BAM::TSequencedBase &data,
							  const TBaseLikelihoods &P_dij_bbar) const noexcept override;
};

struct TWithPMD final: public TModel {
	TBaseLikelihoods P_dij(const BAM::TSequencedBase &, const TBaseLikelihoods &P_dij_bbar) const noexcept override {
		return P_dij_bbar;
	}
	TBaseProbabilities P_bbar(genometools::Base b, const BAM::TSequencedBase &data,
							  const TBaseLikelihoods &P_dij_bbar) const noexcept override;
	TBaseProbabilities P_bbar(genometools::Genotype g, const BAM::TSequencedBase &data,
							  const TBaseLikelihoods &P_dij_bbar) const noexcept override;
};
} // namespace GenotypeLikelihoods::PMD

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
