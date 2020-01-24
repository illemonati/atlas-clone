/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMMODEL_H_
#define TRECALIBRATIONEMMODEL_H_

#include <TRecalibrationEMCovariate.h>
#include "TRecalibrationEMAuxiliaryTools.h"
#include "../TFile.h"

//--------------------------------------------------------------------
// TRecalibrationEMModel
//--------------------------------------------------------------------

class TRecalibrationEMModelNEW{
private:
	size_t _numParameters;
	TRecalibrationEMQualityErrorMap qualErrorMap;

	TLog* logfile;

	//covariates
	TRecalibrationEMCovariateFunction_intercept intercept;
	std::vector<TRecalibrationEMCovariate*> covariates;


	std::vector<TRecalibrationEMCovariateFunction*> activeCovariateFunctions;

	//Newton Raphson Parameters
	//TODO: maybe split into class that can and cannot estimate?
	double _Q, _oldQ;
	arma::mat Jacobian;
	arma::vec F;
	arma::mat JxF;
	unsigned int _numSitesAdded;
	bool _NRconverged;
	bool _NRStepAccepted;

	double _calcEpsilon(const double eta);
	double _calcQ(const int & genotype, TRecalibrationEMReadData & data);

public:
	TRecalibrationEMModelNEW(std::map<std::string, std::string> & modules, TLog* Logfile, bool verbose);
	~TRecalibrationEMModelNEW();

	void createModules(std::map<std::string, std::string> & modules, TRecalibrationEMDataTable* dataTable);
	void createModules(std::map<std::string, std::string> & modules);
	void checkParameterRange(TRecalibrationEMDataTable* dataTable);
	void checkParameterRange(std::vector<uint16_t> & Qualities, uint16_t maxPos);
	int numParameters(){ return _numParameters; };

	void setEMParamsToZero();
	void setQToZero();
	void addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype);
	void addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta);
	double curQ(){ return _Q; };
	bool solveJxF();
	void proposeNewParameters(double & lambda);
	bool acceptProposedParametersBasedOnQ();

	double getSteepestGradient();
	void printJacobianToStdOut();
	void printFToStdOut();
	void printJxFToStdOut();

	double getErrorRate(const TBase & base);
	double getErrorRate(const TRecalibrationEMReadData & data);

};

//--------------------------------------------------------------------
// Global function to create models
//--------------------------------------------------------------------
TRecalibrationEMModel_Base* createTRecalibrationEMModel(std::string modelTag, std::vector<std::string> & values, bool verbose, TLog* logfile);
TRecalibrationEMModel_Base* createTRecalibrationEMModel(std::string modelTag, int maxPos, bool verbose, TLog* logfile);

//--------------------------------------------------------------------
// TRecalibrationEMModels
// Object containing a vector of recal models
//--------------------------------------------------------------------
class TRecalibrationEMModels{
private:
	std::vector<TRecalibrationEMModelNEW*> models;
	unsigned int totNumParameters;
	TRecalibrationEMReadGroupIndex readGroupIndex;
	TLog* logfile;

	void _addModel(std::string & modelTag, std::vector<std::string> & values, bool verbose);
	void _createModelsFromString(std::string & string, TReadGroups & readGroups);
	void _createModelsFromFile(std::string filename, TReadGroups & readGroups, TReadGroupMap & readGroupMap);

	void _writeParameters(TOutputFilePlain & out, const std::string & readGroupName, const int & readGroup, bool isSecondMate);

public:
	TRecalibrationEMModels(int numReadGroups, TLog* Logfile);
	~TRecalibrationEMModels();

	void addSameModelForAllReadGroups(std::string modelTag, std::vector<std::string> & values, bool verbose);
	void addModel(int readGroupId, bool isSecondMate, std::string modelTag, std::vector<std::string> & values, bool verbose);
	void addModel(int readGroupId, bool isSecondMate, std::string modelTag, int maxPos);
	void addModelIfItDoesNotExist(int readGroupId, bool isSecondMate, std::string modelTag, std::vector<int> & Qualities, int maxPos);
	void removeModel(int readGroupId, bool isSecondMate);
	void createModels(std::string string, TReadGroups & readGroups, TReadGroupMap & readGroupMap);

	inline TRecalibrationEMModel_Base* operator[](int index){ return models[index]; };
	int numModels(){ return models.size(); };
	bool modelExists(int readGroupId, bool isSecondMate){ return readGroupIndex.inUse(readGroupId, isSecondMate); };
	inline double calcEpsilon(const TRecalibrationEMReadData & data){
		return models[ readGroupIndex.index(data) ]->getErrorRate(data);
	};

	bool hasReadGroupsWithoutModel();
	void addNoRecalModelIfMissing();
	void reportReadGroupsNotUsed(TReadGroups & readGroups, TReadGroupMap & ReadGroupMap);
	void reportReadGroupsNotUsed(TReadGroups & readGroups);
	void reportReadGroupsConsideredSingleEnd(TReadGroups & readGroups, TReadGroupMap & ReadGroupMap);
	void reportReadGroupsConsideredSingleEnd(TReadGroups & readGroups);
	void warningForMissingReadGroups(TReadGroups & readGroups, TReadGroupMap & ReadGroupMap);
	void warningForMissingReadGroups(TReadGroups & readGroups);

	void setEMParamsToZero();
	void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weight, const double & weightJacobian);
	void setQToZero();
	void addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta);
	void addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype);
	double curQ();
	bool solveJxF();
	void proposeNewParameters(double lambda);
	void scaleParameters();
	unsigned int acceptProposedParametersBasedOnQ();
	void rejectProposedParameters();

	double getSteepestGradient();
	void writeHeader(TOutputFilePlain & out);
	void writeParameters(TOutputFilePlain & out, TReadGroups & readGroups, TReadGroupMap & ReadGroupMap);
	inline double getErrorRate(TBase & base){
		return models[ readGroupIndex.index(base.readGroup, base.isSecondMate) ]->getErrorRate(base);
	};
};


#endif /* TRECALIBRATIONEMMODEL_H_ */
