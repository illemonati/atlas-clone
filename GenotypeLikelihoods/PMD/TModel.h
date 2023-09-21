/*
 * PMD/TModels.h
 */

#ifndef GENOTYPELIKELIHOODS_PMD_MODEL_H_
#define GENOTYPELIKELIHOODS_PMD_MODEL_H_


#include "PMD/TPsi.h"
#include "TAlignment.h"
#include "TGenotypeData.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Types/probability.h"
namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }

namespace GenotypeLikelihoods::PMD {

struct TModel {
	virtual TBaseLikelihoods P_dij(const BAM::TSequencedBase &data,
								   const TBaseLikelihoods &P_dij_bbar) const noexcept          = 0;
	virtual TBaseProbabilities P_bbar(genometools::Base b, const BAM::TSequencedBase &data,
									  const TBaseLikelihoods &P_dij_bbar) const noexcept       = 0;
	virtual TBaseProbabilities P_bbar(genometools::Genotype g, const BAM::TSequencedBase &data,
									  const TBaseLikelihoods &P_dij_bbar) const noexcept       = 0;
	virtual TBaseBaseProbabilities P_b_bbar(genometools::Genotype g, const BAM::TSequencedBase &data,
											const TBaseLikelihoods &P_dij_bbar) const noexcept = 0;
	virtual void simulate(BAM::TAlignment &aln) const                                          = 0;
	virtual TPsi *psi() noexcept                                                               = 0;
	virtual bool hasPMD() const noexcept                                                       = 0;
	virtual BAM::RGInfo::TInfo info() const                                                    = 0;
};

struct TNoPMD final : public TModel {
	TBaseLikelihoods P_dij(const BAM::TSequencedBase &, const TBaseLikelihoods &P_dij_bbar) const noexcept override {
		return P_dij_bbar;
	}
	TBaseProbabilities P_bbar(genometools::Base b, const BAM::TSequencedBase &data,
							  const TBaseLikelihoods &P_dij_bbar) const noexcept override;
	TBaseProbabilities P_bbar(genometools::Genotype g, const BAM::TSequencedBase &data,
							  const TBaseLikelihoods &P_dij_bbar) const noexcept override;
	TBaseBaseProbabilities P_b_bbar(genometools::Genotype g, const BAM::TSequencedBase &data,
									const TBaseLikelihoods &P_dij_bbar) const noexcept override;
	void simulate(BAM::TAlignment &) const override{};
	TPsi *psi() noexcept override {return nullptr;}
	bool hasPMD() const noexcept override {return false;}
	BAM::RGInfo::TInfo info() const override {return {};}
};

class TWithPMD final: public TModel {
	TPsi _psi;
public:
	TWithPMD(std::string_view Psi) : _psi(Psi) {}
	TWithPMD(const BAM::RGInfo::TInfo & info) : _psi(info) {}

	TBaseLikelihoods P_dij(const BAM::TSequencedBase &data, const TBaseLikelihoods &P_dij_bbar) const noexcept override;
	TBaseProbabilities P_bbar(genometools::Base b, const BAM::TSequencedBase &data,
							  const TBaseLikelihoods &P_dij_bbar) const noexcept override;
	TBaseProbabilities P_bbar(genometools::Genotype g, const BAM::TSequencedBase &data,
							  const TBaseLikelihoods &P_dij_bbar) const noexcept override;
	TBaseBaseProbabilities P_b_bbar(genometools::Genotype g, const BAM::TSequencedBase &data,
											const TBaseLikelihoods &P_dij_bbar) const noexcept override;
	void simulate(BAM::TAlignment &aln) const override;
	TPsi *psi() noexcept override {return &_psi;}
	bool hasPMD() const noexcept override {return true;}
	BAM::RGInfo::TInfo info() const override {return _psi.info();}
};
} // namespace GenotypeLikelihoods::PMD

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
