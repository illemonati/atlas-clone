/*
 * SequencingError/TModels.h
 *
 *  Created on: May 14, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_

#include <stddef.h>
#include <stdint.h>
#include <array>
#include <memory>
#include <string>
#include <vector>

#include "TAlignment.h"
#include "coretools/Containers/TStrongArray.h"
#include "genometools/PhredProbabilityTypes.h"
#include "SequencingError/TModel.h"
#include "coretools/Types/probability.h"
#include "TReadGroupInfo.h"

namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------------------
// TModels
// Object containing a vector of TReadGroupModels
//--------------------------------------------------------------------------
using RGModels = coretools::TStrongArray<TModel*, BAM::Mate>;
class TModels {
private:
	std::vector<TWithRecal> _withRecal;
	TNoRecal _noRecal;
	std::vector<RGModels> _pModels;

public:
	void initializeNoRecal(size_t NReadGroups);
	void initialize(std::string_view RecalString, std::string_view RhoString, size_t NReadGroups);
	void initialize(BAM::RGInfo::TReadGroupInfo & RgInfo);

	void pool(const BAM::TReadGroupMap& rgMap);
	void reset(size_t rgID, BAM::Mate mate) noexcept {
		_pModels[rgID][mate] = &_noRecal;
	}

	// access models
	RGModels &RGModel(size_t rgID) noexcept { return _pModels[rgID]; }
	const RGModels &RGModel(size_t rgID) const noexcept { return _pModels[rgID]; }

	TModel& model(const BAM::TSequencedBase &data) noexcept {
		return *RGModel(data.readGroupID)[data.mate()];
	}
	const TModel& model(const BAM::TSequencedBase &data) const noexcept {
		return *RGModel(data.readGroupID)[data.mate()];
	}

	coretools::TView<TWithRecal> unique() noexcept { return _withRecal; }
	coretools::TConstView<TWithRecal> unique() const noexcept { return _withRecal; }

	// operator[]
	const RGModels &operator[](size_t rgID) const noexcept { return RGModel(rgID); }
	RGModels &operator[](size_t rgID) noexcept { return RGModel(rgID); }
	const TModel &operator[](const BAM::TSequencedBase &data) const noexcept { return model(data); }
	TModel &operator[](const BAM::TSequencedBase &data) noexcept { return model(data); }

	bool recalibrates() const noexcept { return !_withRecal.empty(); }

	// calculate error rates
	coretools::Probability errorRate(const BAM::TSequencedBase &data) const noexcept {
		return model(data).errorRate(data);
	}
	genometools::PhredIntProbability phredInt(const BAM::TSequencedBase &data) const noexcept {
		return model(data).phredInt(data);
	}
	void recalibrate(BAM::TSequencedBase &data) const noexcept { data.recalibratedQualityAsPhredInt = phredInt(data); }
	void recalibrate(std::vector<BAM::TSequencedBase> &datas) const noexcept;
	TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &data) const noexcept {
		return model(data).baseLikelihoods(data);
	}

	void addToRGInfo(BAM::RGInfo::TReadGroupInfo & RgInfo) const;
};
} // namespace SequencingError
}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
