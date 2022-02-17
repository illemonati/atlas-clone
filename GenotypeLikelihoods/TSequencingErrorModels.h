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

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------------------
// TSequencingErrorModelsOneReadGroup
//--------------------------------------------------------------------------
class TModelEntry {
private:
	std::shared_ptr<TModelRecal> _recalModel;
	static TModelNoRecal _noRecalModel;

public:
	TModelEntry()  = default;
	~TModelEntry() = default;

	void addModel(const TModelDefinition &ModelDef);
	const TModel &model() const;
	TModelRecal *getPointerToRecalModel();
	std::shared_ptr<TModelRecal> &getSharedPointerToRecalModel();
	bool recalibrates() const { return _recalModel->recalibrates(); };

	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const;
	void fillBaseLikelihoods(const BAM::TSequencedBase &base, TBaseLikelihoods &baseLikelihoods) const;
};

class TModelsOneReadGroup {
private:
	std::array<TModelEntry, 2> _models;

public:
	TModelsOneReadGroup()  = default;
	~TModelsOneReadGroup() = default;

	void addRecalModel(const TModelDefinition &ModelDef, const bool &IsSecondMate);

	const TModel &operator[](const bool &IsSecondMate) const;
	TModelRecal *getPointerToRecalModel(const bool &IsSecondMate);
	std::shared_ptr<TModelRecal> &getSharedPointerToRecalModel(const bool &IsSecondMate);
	bool recalibrates() const;

	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const;
	void fillBaseLikelihoods(const BAM::TSequencedBase &base, TBaseLikelihoods &baseLikelihoods) const;
};

//--------------------------------------------------------------------------
// TSequencingErrorModels
// Object containing a vector of TSequencingErrorModelsOneReadGroup
//--------------------------------------------------------------------------
class TModels {
private:
	// models
	std::vector<std::array<TModelEntry, 2>> _models;

	void _addModel(TModelDefinition &modelDef);
	void _addNoRecalModelIfMissing();

public:
	bool recalStringIsLikelyAModel(const std::string &RecalString);
	void initialize(const std::string &RecalString, const std::string &RhoString, const BAM::TReadGroups &ReadGroups,
			coretools::TLog *Logfile);
	void initializeFromFile(const std::string &Filename, const BAM::TReadGroups &ReadGroups, coretools::TLog *Logfile);
	void checkReadGroups(const BAM::TReadGroups &ReadGroups, std::vector<uint16_t> &ReadGroupsWithoutRecal,
			     std::vector<uint16_t> &ReadGroupsLikelySingleEnd);

	// access models
	TModelEntry &operator()(uint16_t ReadGroupIndex, bool IsSecondMate) noexcept { return _models[ReadGroupIndex][IsSecondMate]; };
	size_t numModels() { return _models.size(); };
	bool recalibrationChangesQualities() const noexcept { return !_models.empty(); }
	bool recalibrates(uint16_t ReadGroupIndex) const noexcept {
		return _models.size() > ReadGroupIndex &&
		       (_models[ReadGroupIndex][0].recalibrates() || _models[ReadGroupIndex][1].recalibrates());
	}

	// calculate error rates
	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const;
	void recalibrate(BAM::TSequencedBase &base) const;
	void recalibrate(std::vector<BAM::TSequencedBase> &bases, const uint16_t length) const; // TODO: remove
	void recalibrate(std::vector<BAM::TSequencedBase> &bases) const;
	void fillBaseLikelihoods(const BAM::TSequencedBase &base, TBaseLikelihoods &baseLikelihoods) const;

	void writeRecalFile(const BAM::TReadGroups &ReadGroups, const std::string & Filename) const;
	void print() const;
};
} // namespace SequencingError
}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
