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

//--------------------------------------------------------------------
// TSequencingErrorModels
// Object containing a vector of TSequencingErrorModel, generally one per read group and mate
//--------------------------------------------------------------------
class TSequencingErrorModels{
private:
	TLog* logfile;
	BAM::TReadGroups* readGroups;
	TSequencingErrorRho defaultRho;

	//models
	bool doRecalibration;
	BAM::TReadGroupMap* readGroupMap;
	std::vector<TSequencingErrorModel> models;
	TRecalibrationEMReadGroupIndex readGroupIndex;
	unsigned int totNumParameters;

	void _init(BAM::TReadGroups* ReadGroups,  BAM::TReadGroupMap* ReadGroupMap, TLog* Logfile);

	void _readRecalFile(const std::string & filename, std::vector<TSequencingErrorModelDefinition> & modelDefs);
	void _addModel(TSequencingErrorModelDefinition & modelDef);
	void _addNoRecalModelIfMissing();

	void _writeParameters(TOutputFile & out, const std::string & readGroupName, const int & readGroup, bool isSecondMate) const;

public:
	BAM::TQualityMap qualMap; //make available to others
	TSequencingErrorModels();

	//add model for recalibration: no dataTable provided
	void createModels(const std::string & filename, BAM::TReadGroups* ReadGroups,  BAM::TReadGroupMap* ReadGroupMap, TLog* Logfile);
	void createEmptyModels(BAM::TReadGroups* ReadGroups,  BAM::TReadGroupMap* ReadGroupMap, TLog* Logfile);

	//add models for estimation: dataTable provided
	void prepareModelsForEstimation(BAM::TReadGroups* ReadGroups,  BAM::TReadGroupMap* ReadGroupMap, TLog* Logfile);
	void addModel(TSequencingErrorModelDefinition & covariates, TRecalibrationEMDataTable* dataTable);
	void addModelsFromFile(std::string filename, TRecalibrationEMDataTables* dataTables);
	void removeModel(int readGroupId, bool isSecondMate);

	int numModels(){ return models.size(); };
	bool modelExists(uint16_t readGroupId, bool isSecondMate){ return readGroupIndex.inUse(readGroupId, isSecondMate); };
	bool modelExists(TSequencingErrorModelDefinition & def){ return readGroupIndex.inUse(def.readGroupId, def.isSecondMate); };
	bool recalibrationChangesQualities() const{ return doRecalibration; };

	bool hasReadGroupsWithoutModel() const;
	void reportReadGroupsNotUsed() const;
	void reportReadGroupsConsideredSingleEnd() const;
	void warningForMissingReadGroups() const;

	//calculate error rates
	double getErrorRate(const BAM::TBase & base) const;
	uint8_t getPhredInt(const BAM::TBase & base) const;
	void recalibrate(BAM::TBase & base) const;
	void recalibrate(std::vector<BAM::TBase> & bases, const uint16_t  length) const; //TODO: remove
	void recalibrate(std::vector<BAM::TBase> & bases) const;
	void calculateBaseLikelihoods(const BAM::TBase & base, TBaseData & baseLikelihoods) const;

	//functions to estimate rho
	void prepareRhoEstimationFromEMWeights();
	void addBaseForRhoEstimation(BAM::TBase & base, const TBaseData & EMWeights);
	void estimateRho();

	//functions to estimate beta
	void setNewtonRaphsonParamsToZero();
	void addToFandJacobian(const BAM::TBase & base, const TBaseData & EM_weights_bbar_given_d);
	void setQToZero();
	void addToQ(const BAM::TBase & base, const TBaseData & EM_weights_bbar_given_d);
	double curQ();
	bool solveJxF();
	void proposeNewParameters(double lambda);
	void scaleParameters();
	unsigned int acceptProposedParametersBasedOnQ();
	void adjustParametersPostEstimation();
	double getSteepestGradient();

	void writeRecalFile(const std::string filename) const;
	void print() const;
};

}; //end namespace


#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
