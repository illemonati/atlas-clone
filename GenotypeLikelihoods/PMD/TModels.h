/*
 * PMD/TModels.h
 *
 */

#ifndef GENOTYPELIKELIHOODS_PMD_MODELS_H_
#define GENOTYPELIKELIHOODS_PMD_MODELS_H_

#include <vector>

#include "TReadGroupMap.h"
#include "genometools/GenotypeTypes.h"

#include "TGenotypeData.h"
#include "TSequencedBase.h"
#include "TModel.h"

namespace BAM::RGInfo {class TReadGroupInfo;}

namespace GenotypeLikelihoods::PMD {

class TModels {
private:
	std::vector<TWithPMD> _withPMD;
	TNoPMD _noPMD;
	std::vector<TModel*> _pModels;

public:
	void initialize(size_t NReadGroups, std::string_view PMDString = "");
	void initialize(const BAM::RGInfo::TReadGroupInfo & RgInfo);

	void pool(const BAM::TReadGroupMap& rgMap);
	void reset(size_t rgID) noexcept {
		_pModels[rgID] = &_noPMD;
	}

	// access models
	TModel &model(size_t rgID) noexcept { return *_pModels[rgID]; }
	const TModel &model(size_t rgID) const noexcept { return *_pModels[rgID]; }

	TModel& model(const BAM::TSequencedBase &data) noexcept {
		return model(data.readGroupID);
	}
	const TModel& model(const BAM::TSequencedBase &data) const noexcept {
		return model(data.readGroupID);
	}

	bool hasPMD() const noexcept { return !_withPMD.empty(); }

	TBaseLikelihoods P_dij(const BAM::TSequencedBase &data, const TBaseLikelihoods &P_dij_bbar) const noexcept {
		return model(data).P_dij(data, P_dij_bbar);
	}

	TBaseProbabilities P_bbar(genometools::Base b, const BAM::TSequencedBase &data,
							  const TBaseLikelihoods &P_dij_bbar) const noexcept {
		return model(data).P_bbar(b, data, P_dij_bbar);
	}

	TBaseProbabilities P_bbar(genometools::Genotype g, const BAM::TSequencedBase &data,
							  const TBaseLikelihoods &P_dij_bbar) const noexcept {
		return model(data).P_bbar(g, data, P_dij_bbar);
	}

	TBaseBaseProbabilities P_b_bar(genometools::Genotype g, const BAM::TSequencedBase &data,
										   const TBaseLikelihoods &P_dij_bbar) const noexcept {
		return model(data).P_b_bbar(g, data, P_dij_bbar);
	}

	void addToRGInfo(BAM::RGInfo::TReadGroupInfo & RgInfo) const;
};
} // namespace GenotypeLikelihoods::PMD

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
