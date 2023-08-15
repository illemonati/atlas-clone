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
#include "genometools/PhredProbabilityTypes.h"
#include "SequencingError/TModel.h"
#include "coretools/Types/probability.h"
#include "TReadGroupInfo.h"

namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }

namespace GenotypeLikelihoods {
namespace SequencingError {

//-------------------------------------
// TReadGroupModels
// Object containing a TModel for the first and second mate
//-------------------------------------
class TReadGroupModels{
private:
	std::array<std::unique_ptr<TModel>, 2> _models;

public:
	TReadGroupModels();
	TReadGroupModels(std::string_view RecalString, std::string_view RhoString);
	TReadGroupModels(std::string_view RecalString1, std::string_view RhoString1, std::string_view RecalString2, std::string_view RhoString2);
	TReadGroupModels(const BAM::RGInfo::TReadGroupInfoEntry & Info);

	void initialize(size_t mate, std::string_view RecalString, std::string_view RhoString);

	TModel& operator[](bool isSecondMate) noexcept {
		return *_models[isSecondMate].get();
	}

	const TModel& operator[](bool isSecondMate) const noexcept {
		return *_models[isSecondMate].get();
	}

	bool recalibrates() const  noexcept {
		if (_models[0]->recalibrates() || _models[1]->recalibrates()) return true;
		return false;
	}

	void simulate(BAM::TAlignment & Alignment) const;

	BAM::RGInfo::TInfo info() const;
};

//--------------------------------------------------------------------------
// TModels
// Object containing a vector of TReadGroupModels
//--------------------------------------------------------------------------
class TModels {
private:
	std::vector<TReadGroupModels> _models;
	std::vector<TReadGroupModels*> _pModels;
	bool _recalibrates = false;

	TReadGroupModels &_model(size_t rgID) noexcept { return *_pModels[rgID]; }
	const TReadGroupModels &_model(size_t rgID) const noexcept { return *_pModels[rgID]; }

public:
	void initializeNoRecal(const BAM::TReadGroups &ReadGroups);
	void initialize(std::string_view RecalString, std::string_view RhoString, const BAM::TReadGroups &ReadGroups);
	void initialize(const std::vector<std::string> & RecalStringPerReadGroup, const std::vector<std::string> & RhoStringPerReadGroup, const BAM::TReadGroups &ReadGroups);
	void initializeFromFile(std::string_view Filename, const BAM::TReadGroups &ReadGroups);
	void initialize(BAM::RGInfo::TReadGroupInfo & RgInfo);
	void checkReadGroups(const BAM::TReadGroups &ReadGroups, std::vector<size_t> &ReadGroupsWithoutRecal,
			     std::vector<size_t> &ReadGroupsLikelySingleEnd) const noexcept;

	std::vector<TReadGroupModels> forget();
	void remember(std::vector<TReadGroupModels>& forgottenModels);

	void pool(const BAM::TReadGroupMap& rgMap);
	void unPool();

	// access models
	TModel& model(const BAM::TSequencedBase &data) noexcept {
		return _model(data.readGroupID)[data.isSecondMate()];
	}
	const TModel& model(const BAM::TSequencedBase &data) const noexcept {
		return _model(data.readGroupID)[data.isSecondMate()];
	}

	const TReadGroupModels &operator[](size_t rgID) const noexcept { return _model(rgID); }
	TReadGroupModels &operator[](size_t rgID) noexcept { return _model(rgID); }
	const TModel &operator[](const BAM::TSequencedBase &data) const noexcept { return model(data); }
	TModel &operator[](const BAM::TSequencedBase &data) noexcept { return model(data); }

	bool recalibrates() const noexcept { return _recalibrates; }

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

	void writeRecalFile(const BAM::TReadGroups &ReadGroups, std::string_view Filename) const;

	void addToRGInfo(BAM::RGInfo::TReadGroupInfo & RgInfo) const;

	std::string getRhoDefinition() {
		if (_models.empty()) return std::string{};
		return "TODO";
	}

	std::string getModelsDefinition() {
		if (_models.empty()) return std::string{};
		return "TODO";
	}
};
} // namespace SequencingError
}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
