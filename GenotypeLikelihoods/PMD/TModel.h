/*
 * PMD/TModels.h
 */

#ifndef GENOTYPELIKELIHOODS_PMD_MODEL_H_
#define GENOTYPELIKELIHOODS_PMD_MODEL_H_

#include "PMD/TPsi.h"
#include "TAlignment.h"
#include "coretools/Types/probability.h"
#include "genometools/Genotypes/Containers.h"

namespace BAM {
class TReadGroups;
struct TSequencedData;
} // namespace BAM

namespace GenotypeLikelihoods::PMD {

struct TModel {
	virtual genometools::TBaseLikelihoods P_dij(const BAM::TSequencedData &data,
												const genometools::TBaseLikelihoods &P_dij_bbar) const noexcept    = 0;
	virtual genometools::TBaseProbabilities P_bbar(genometools::Base b, const BAM::TSequencedData &data,
												   const genometools::TBaseLikelihoods &P_dij_bbar) const noexcept = 0;
	virtual genometools::TBaseProbabilities P_bbar(genometools::Genotype g, const BAM::TSequencedData &data,
												   const genometools::TBaseLikelihoods &P_dij_bbar) const noexcept = 0;
	virtual TBaseBaseProbabilities P_b_bbar(genometools::Genotype g, const BAM::TSequencedData &data,
											const genometools::TBaseLikelihoods &P_dij_bbar) const noexcept        = 0;
	virtual void simulate(BAM::TAlignment &aln) const                                                              = 0;
	virtual TPsi *psi() noexcept                                                                                   = 0;
	virtual bool hasPMD() const noexcept                                                                           = 0;
	virtual BAM::RGInfo::TInfo info() const                                                                        = 0;
	virtual coretools::Probability prob(Type T, const BAM::TSequencedData &Data) const                             = 0;
};

struct TNoPMD final : public TModel {
	genometools::TBaseLikelihoods P_dij(const BAM::TSequencedData &, const genometools::TBaseLikelihoods &P_dij_bbar) const noexcept override {
		return P_dij_bbar;
	}
	genometools::TBaseProbabilities P_bbar(genometools::Base b, const BAM::TSequencedData &data,
							  const genometools::TBaseLikelihoods &P_dij_bbar) const noexcept override;
	genometools::TBaseProbabilities P_bbar(genometools::Genotype g, const BAM::TSequencedData &data,
							  const genometools::TBaseLikelihoods &P_dij_bbar) const noexcept override;
	TBaseBaseProbabilities P_b_bbar(genometools::Genotype g, const BAM::TSequencedData &data,
									const genometools::TBaseLikelihoods &P_dij_bbar) const noexcept override;
	void simulate(BAM::TAlignment &) const override{};
	TPsi *psi() noexcept override {return nullptr;}
	bool hasPMD() const noexcept override {return false;}
	BAM::RGInfo::TInfo info() const override { return "-"; }
	coretools::Probability prob(Type, const BAM::TSequencedData &) const override { return coretools::Probability{0.}; }
};

class TWithPMD final: public TModel {
	TPsi _psi;
public:
	TWithPMD() = default;
	TWithPMD(const BAM::RGInfo::TInfo & info) : _psi(info) {}

	genometools::TBaseLikelihoods P_dij(const BAM::TSequencedData &data, const genometools::TBaseLikelihoods &P_dij_bbar) const noexcept override;
	genometools::TBaseProbabilities P_bbar(genometools::Base b, const BAM::TSequencedData &data,
							  const genometools::TBaseLikelihoods &P_dij_bbar) const noexcept override;
	genometools::TBaseProbabilities P_bbar(genometools::Genotype g, const BAM::TSequencedData &data,
							  const genometools::TBaseLikelihoods &P_dij_bbar) const noexcept override;
	TBaseBaseProbabilities P_b_bbar(genometools::Genotype g, const BAM::TSequencedData &data,
											const genometools::TBaseLikelihoods &P_dij_bbar) const noexcept override;
	void simulate(BAM::TAlignment &aln) const override;
	TPsi *psi() noexcept override {return &_psi;}
	bool hasPMD() const noexcept override {return true;}
	BAM::RGInfo::TInfo info() const override {return _psi.info();}
	coretools::Probability prob(Type T, const BAM::TSequencedData &Data) const override {return _psi.prob(T, Data);}
};
} // namespace GenotypeLikelihoods::PMD

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
