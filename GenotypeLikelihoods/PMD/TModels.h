/*
 * PMD/TModels.h
 *
 */

#ifndef GENOTYPELIKELIHOODS_PMD_MODELS_H_
#define GENOTYPELIKELIHOODS_PMD_MODELS_H_

#include "TGenotypeData.h"
#include "TReadGroupInfo.h"
#include "TSequencedBase.h"
#include <string>
#include <vector>

namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }

namespace GenotypeLikelihoods::PMD {

struct TModel {
	TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &) const noexcept {
		return TBaseLikelihoods{};
	}
};
struct TWithPMD final: public TModel {};
struct TNoPMD final: public TModel {};

class TModels {
private:
	std::vector<TWithPMD> _withPMD;
	TNoPMD _noPMD;
	std::vector<TModel*> _pModels;

public:
	void initialize(size_t NReadGroups, std::string_view PMDString = "");
	void initialize(BAM::RGInfo::TReadGroupInfo & RgInfo);

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

	bool recalibrates() const noexcept { return !_withPMD.empty(); }

	TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &data) const noexcept {
		return model(data).baseLikelihoods(data);
	}

	void addToRGInfo(BAM::RGInfo::TReadGroupInfo & RgInfo) const;
};
} // namespace GenotypeLikelihoods::PMD

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
