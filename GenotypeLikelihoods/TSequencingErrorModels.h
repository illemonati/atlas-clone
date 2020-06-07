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

	void _readRecalFile(const std::string filename, std::vector<TSequencingErrorModelDefinition> & modelDefs);

	void _createModelsFromString(const std::string & s);
	void _createModelsFromFile(std::string filename);
	void _addModel(const uint16_t readGroupId, const bool isSecondMate, TSequencingErrorCovariateDefinition & covariateFunctions);
	void _addNoRecalModelIfMissing();

	void _writeParameters(TOutputFile & out, const std::string & readGroupName, const int & readGroup, bool isSecondMate) const;

public:
	TQualityMap qualMap; //make available to others
	TSequencingErrorModels();

	//add model for recalibration: no dataTable provided
	void createModels(std::string string, BAM::TReadGroups* ReadGroups,  BAM::TReadGroupMap* ReadGroupMap, TLog* Logfile);
	void createEmptyModels(BAM::TReadGroups* ReadGroups,  BAM::TReadGroupMap* ReadGroupMap, TLog* Logfile);

	//add models for estimation: dataTable provided
	void prepareModelsForEstimation(BAM::TReadGroups* ReadGroups,  BAM::TReadGroupMap* ReadGroupMap, TLog* Logfile);
	void addModel(const uint16_t readGroupId, const bool isSecondMate, TSequencingErrorCovariateDefinition & covariates, TRecalibrationEMDataTable* dataTable);
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
	double getErrorRate(const TRecalibrationEMReadData & data) const;
	double getErrorRate(const TBase & base) const;
	uint8_t getPhredInt(const TBase & base) const;
	void recalibrate(TBase & base) const;
	void recalibrate(std::vector<TBase> & bases, const uint16_t  length) const; //TODO: remove
	void recalibrate(std::vector<TBase> & bases) const;
	void calculateBaseLikelihoods(const TBase & base, TBaseData & baseLikelihoods) const;

	//function to estimate
	void setEMParamsToZero();
	void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	void setQToZero();
	void addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta);
	void addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype);
	double curQ();
	bool solveJxF();
	void proposeNewParameters(double lambda);
	void scaleParameters();
	unsigned int acceptProposedParametersBasedOnQ();
	void adjustParametersPostEstimation();
	double getSteepestGradient();

	void writeRecalFile(const std::string filename) const;
};

}; //end namespace


#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
