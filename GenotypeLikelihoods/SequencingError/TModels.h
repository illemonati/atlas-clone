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

#include "PhredProbabilityTypes.h"
#include "SequencingError/TModel.h"
#include "probability.h"

namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------------------
// TModel
// Object containing a vector of TModel
//--------------------------------------------------------------------------
class TModels {
private:
	std::vector<std::array<std::unique_ptr<TModel>, 2>> _models;
public:
	void initializeNoRecal(const BAM::TReadGroups &ReadGroups);
	void initialize(const std::string &RecalString, const std::string &RhoString, const BAM::TReadGroups &ReadGroups);
	void initializeFromFile(const std::string &Filename, const BAM::TReadGroups &ReadGroups);
	void checkReadGroups(const BAM::TReadGroups &ReadGroups, std::vector<uint16_t> &ReadGroupsWithoutRecal,
			     std::vector<uint16_t> &ReadGroupsLikelySingleEnd) const noexcept;

	std::vector<std::array<std::unique_ptr<TModel>, 2>> forget();
	void remember(std::vector<std::array<std::unique_ptr<TModel>, 2>>& forgottenModels);

	// access models
	TModel &operator()(uint16_t ReadGroupIndex, bool IsSecondMate) noexcept {
		return *_models[ReadGroupIndex][IsSecondMate];
	}

	const TModel &operator()(uint16_t ReadGroupIndex, bool IsSecondMate) const noexcept {
		return *_models[ReadGroupIndex][IsSecondMate];
	}

	bool recalibrationChangesQualities() const noexcept {
		for (const auto& m: _models) {
			if (m[0]->recalibrates() || m[1]->recalibrates()) return true;
		}
		return false;
	}

	// calculate error rates
	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const noexcept;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const noexcept;
	void recalibrate(BAM::TSequencedBase &base) const noexcept;
	void recalibrate(std::vector<BAM::TSequencedBase> &bases) const noexcept;
	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept;

	void writeRecalFile(const BAM::TReadGroups &ReadGroups, const std::string & Filename) const;
};
} // namespace SequencingError
}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
