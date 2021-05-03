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
class TSequencingErrorModelsOneReadGroup{
private:
	//stores pointers to recal models for estimation
	std::shared_ptr<TSequencingErrorModelRecal> _recalModelFirstMate;
	std::shared_ptr<TSequencingErrorModelRecal> _recalModelSecondMate;

	//and pointers to general models for ease of use
	std::shared_ptr<TSequencingErrorModel> _modelFirstMate;
	std::shared_ptr<TSequencingErrorModel> _modelSecondMate;

public:
	TSequencingErrorModelsOneReadGroup(std::shared_ptr<TSequencingErrorModelNoRecal> & NoRecalModel);
	~TSequencingErrorModelsOneReadGroup() = default;

	void addRecalModel(const TSequencingErrorModelDefinition & ModelDef, const bool & IsSecondMate);

	const TSequencingErrorModel& operator[](const bool & IsSecondMate) const;
	const std::shared_ptr<TSequencingErrorModelRecal>& getPointerToRecalModel(bool IsSecondMate) const;

	double getErrorRate(const BAM::TBase & base, const BAM::TQualityMap & qualMap) const;
	uint8_t getPhredInt(const BAM::TBase & base, const BAM::TQualityMap & qualMap) const;
	void fillBaseLikelihoods(const BAM::TBase & base, const BAM::TQualityMap & qualMap, TBaseData & baseLikelihoods) const;
};

//--------------------------------------------------------------------------
// TSequencingErrorModels
// Object containing a vector of TSequencingErrorModelsOneReadGroup
//--------------------------------------------------------------------------
class TSequencingErrorModels{
private:
	BAM::TQualityMap qualMap;

	//models
	bool doRecalibration;
	std::shared_ptr<TSequencingErrorModelNoRecal> _noRecalModel;
	std::vector<TSequencingErrorModelsOneReadGroup> models;

	void _initializeFromFile(const std::string & filename, std::vector<TSequencingErrorModelDefinition> & modelDefs);
	void _addModel(TSequencingErrorModelDefinition & modelDef);
	void _addNoRecalModelIfMissing();

	void _writeParameters(TOutputFile & out, const std::string & readGroupName, const int & readGroup, bool isSecondMate) const;

public:

	TSequencingErrorModels();

	void initialize(const std::string & RecalString, const std::string & RhoString, const BAM::TReadGroups & ReadGroups, TLog* Logfile);
	void initializeFromFile(const std::string & Filename, const BAM::TReadGroups & ReadGroups, TLog* Logfile, std::vector<uint16_t> & ReadGroupsWithoutRecal, std::vector<uint16_t> & ReadGroupsLikelySingleEnd);

	//access models
	const TSequencingErrorModelsOneReadGroup& operator[](const uint16_t & ReadGroupIndex){ return models[ReadGroupIndex]; };

	int numModels(){ return models.size(); };
	bool recalibrationChangesQualities() const{ return doRecalibration; };

	//calculate error rates
	double getErrorRate(const BAM::TBase & base) const;
	uint8_t getPhredInt(const BAM::TBase & base) const;
	void recalibrate(BAM::TBase & base) const;
	void recalibrate(std::vector<BAM::TBase> & bases, const uint16_t  length) const; //TODO: remove
	void recalibrate(std::vector<BAM::TBase> & bases) const;
	void fillBaseLikelihoods(const BAM::TBase & base, TBaseData & baseLikelihoods) const;

	void writeRecalFile(const std::string filename) const;
	void print() const;
};

}; //end namespace


#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
