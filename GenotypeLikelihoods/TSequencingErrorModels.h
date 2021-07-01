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

#include "TSequencingErrorModel.h"
#include "TFile.h"
#include "auxiliaryTools.h"

namespace GenotypeLikelihoods{


//--------------------------------------------------------------------------
// TSequencingErrorModelsOneReadGroup
//--------------------------------------------------------------------------
class TSequencingErrorModelEntry{
private:
	std::shared_ptr<TSequencingErrorModelRecal> _recalModel;
	static TSequencingErrorModelNoRecal _noRecalModel;

public:
	TSequencingErrorModelEntry() = default;
	~TSequencingErrorModelEntry() = default;

	void addModel(const TSequencingErrorModelDefinition & ModelDef);
	const TSequencingErrorModel& model() const;
	TSequencingErrorModelRecal* getPointerToRecalModel();
	std::shared_ptr<TSequencingErrorModelRecal>& getSharedPointerToRecalModel();

	Probability getErrorRate(const BAM::TSequencedBase & base) const;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase & base) const;
	void fillBaseLikelihoods(const BAM::TSequencedBase & base, TBaseLikelihoods & baseLikelihoods) const;
};

class TSequencingErrorModelsOneReadGroup{
private:
	std::array<TSequencingErrorModelEntry,2> _models;

public:
	TSequencingErrorModelsOneReadGroup();
	~TSequencingErrorModelsOneReadGroup() = default;

	void addRecalModel(const TSequencingErrorModelDefinition & ModelDef, const bool & IsSecondMate);

	const TSequencingErrorModel& operator[](const bool & IsSecondMate) const;
	TSequencingErrorModelRecal* getPointerToRecalModel(const bool & IsSecondMate);
	std::shared_ptr<TSequencingErrorModelRecal>& getSharedPointerToRecalModel(const bool & IsSecondMate);

	Probability getErrorRate(const BAM::TSequencedBase & base) const;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase & base) const;
	void fillBaseLikelihoods(const BAM::TSequencedBase & base, TBaseLikelihoods & baseLikelihoods) const;
};

//--------------------------------------------------------------------------
// TSequencingErrorModels
// Object containing a vector of TSequencingErrorModelsOneReadGroup
//--------------------------------------------------------------------------
class TSequencingErrorModels{
private:
	//models
	bool _doRecalibration;
	TSequencingErrorModelNoRecal _noRecalModel;
	std::vector<TSequencingErrorModelsOneReadGroup> _models;

	void _addModel(TSequencingErrorModelDefinition & modelDef);
	void _addNoRecalModelIfMissing();

public:
	TSequencingErrorModels();

	bool recalStringIsLikelyAModel(const std::string & RecalString);
	void initialize(const std::string & RecalString, const std::string & RhoString, const BAM::TReadGroups & ReadGroups, coretools::TLog* Logfile);
	void initializeFromFile(const std::string & Filename, const BAM::TReadGroups & ReadGroups, coretools::TLog* Logfile);
	void checkReadGroups(const BAM::TReadGroups & ReadGroups, std::vector<uint16_t> & ReadGroupsWithoutRecal, std::vector<uint16_t> & ReadGroupsLikelySingleEnd);

	//access models
	TSequencingErrorModelsOneReadGroup& operator[](const uint16_t & ReadGroupIndex){ return _models[ReadGroupIndex]; };
	size_t numModels(){ return _models.size(); };
	bool recalibrationChangesQualities() const{ return _doRecalibration; };

	//calculate error rates
	Probability getErrorRate(const BAM::TSequencedBase & base) const;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase & base) const;
	void recalibrate(BAM::TSequencedBase & base) const;
	void recalibrate(std::vector<BAM::TSequencedBase> & bases, const uint16_t length) const; //TODO: remove
	void recalibrate(std::vector<BAM::TSequencedBase> & bases) const;
	void fillBaseLikelihoods(const BAM::TSequencedBase & base, TBaseLikelihoods & baseLikelihoods) const;

	void writeRecalFile(const BAM::TReadGroups & ReadGroups, const std::string Filename) const;
	void print() const;
};

}; //end namespace


#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
