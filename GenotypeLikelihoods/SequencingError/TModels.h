/*
 * SequencingError/TModels.h
 *
 *  Created on: May 14, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_

#include <string>
#include <vector>

#include "SequencingError/TModel.h"
#include "TAlignment.h"
#include "TReadGroupInfo.h"
#include "coretools/Containers/TStrongArray.h"
#include "genometools/PhredProbabilityTypes.h"

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
	void _initializeNoRecal(size_t NReadGroups);

public:
	void initialize(size_t NReadGroups, std::string_view RecalString = "", std::string_view RhoString = "");
	void initialize(const BAM::RGInfo::TReadGroupInfo & RgInfo);

	void pool(const BAM::TReadGroupMap& rgMap);
	void reset(size_t rgID, BAM::Mate mate) noexcept { _pModels[rgID][mate] = &_noRecal; }

	// access models
	RGModels &RGModel(size_t rgID) noexcept { return _pModels[rgID]; }
	const RGModels &RGModel(size_t rgID) const noexcept { return _pModels[rgID]; }

	TModel& model(const BAM::TSequencedBase &data) noexcept {
		return *RGModel(data.readGroupID)[data.mate()];
	}
	const TModel& model(const BAM::TSequencedBase &data) const noexcept {
		return *RGModel(data.readGroupID)[data.mate()];
	}

	bool recalibrates() const noexcept { return !_withRecal.empty(); }

	// calculate error rates
	genometools::PhredIntProbability phredInt(const BAM::TSequencedBase &data) const noexcept {
		return model(data).phredInt(data);
	}
	void recalibrate(BAM::TAlignment &aln) const noexcept { RGModel(aln.readGroupId())[aln.mate()]->recalibrate(aln); };
	TBaseLikelihoods P_dij(const BAM::TSequencedBase &data) const noexcept {
		return model(data).P_dij(data);
	}

	void addToRGInfo(BAM::RGInfo::TReadGroupInfo & RgInfo) const;
};
} // namespace SequencingError
}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
