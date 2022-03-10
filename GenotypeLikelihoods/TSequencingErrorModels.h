/*
 * TSequencingErrorModels.h
 *
 *  Created on: May 14, 2020
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_

/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TFile.h"
#include "TSequencingErrorModel.h"
#include "auxiliaryTools.h"
#include <memory>

namespace GenotypeLikelihoods {
namespace SequencingError {

class TModelsMates {};

//--------------------------------------------------------------------------
// TModel
// Object containing a vector of TModel
//--------------------------------------------------------------------------
class TModels {
private:
	std::vector<std::array<std::shared_ptr<TModelRecal>, 2>> _models;
	TModelNoRecal _noRecal;
public:
	bool recalStringIsLikelyAModel(const std::string &RecalString) const noexcept;
	void initialize(const std::string &RecalString, const std::string &RhoString, const BAM::TReadGroups &ReadGroups);
	void initializeFromFile(const std::string &Filename, const BAM::TReadGroups &ReadGroups);
	void checkReadGroups(const BAM::TReadGroups &ReadGroups, std::vector<uint16_t> &ReadGroupsWithoutRecal,
			     std::vector<uint16_t> &ReadGroupsLikelySingleEnd) const noexcept;

	// access models
	TModel &operator()(uint16_t ReadGroupIndex, bool IsSecondMate) noexcept {
		if (_models.size() < ReadGroupIndex && _models[ReadGroupIndex][IsSecondMate])
			return *_models[ReadGroupIndex][IsSecondMate];
		return _noRecal;
	}

	const TModel &operator()(uint16_t ReadGroupIndex, bool IsSecondMate) const noexcept {
		if (_models.size() < ReadGroupIndex && _models[ReadGroupIndex][IsSecondMate])
			return *_models[ReadGroupIndex][IsSecondMate];
		return _noRecal;
	}

	std::shared_ptr<TModelRecal> getRecal(uint16_t ReadGroupIndex, bool IsSecondMate) noexcept {
		return _models[ReadGroupIndex][IsSecondMate];
	}

	size_t numModels() const noexcept { return _models.size(); };
	bool recalibrationChangesQualities() const noexcept { return !_models.empty(); }

	// calculate error rates
	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const noexcept;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const noexcept;
	void recalibrate(BAM::TSequencedBase &base) const noexcept;
	void recalibrate(std::vector<BAM::TSequencedBase> &bases) const noexcept;
	void fillBaseLikelihoods(const BAM::TSequencedBase &base, TBaseLikelihoods &baseLikelihoods) const noexcept;

	void writeRecalFile(const BAM::TReadGroups &ReadGroups, const std::string & Filename) const;
	void print() const;
};
} // namespace SequencingError
}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
