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
// TRecalibrationEMModelCovariateDefinition
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
class TRecalibrationEMModelCovariateDefinition{
private:
	std::map<std::string, std::string> covariateFunctions;  //<covariate, function>

public:
	TRecalibrationEMModelCovariateDefinition(){};
	TRecalibrationEMModelCovariateDefinition(const std::string modelString, std::string & error){
		parse(modelString, error);
	};

	bool parse(const std::string & modelString, std::string & error);
	void add(const std::string covariate, const std::string function){
		covariateFunctions.insert(std::pair<std::string, std::string>(covariate, function));
	};
	size_t size(){ return covariateFunctions.size(); };
	std::map<std::string, std::string>::iterator begin(){ return covariateFunctions.begin(); };
	std::map<std::string, std::string>::iterator end(){ return covariateFunctions.end(); };
	std::string getModelString();
};

//--------------------------------------------------------------------
// TRecalibrationEMModelDefinition
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
class TRecalibrationEMModelDefinition{
public:
	std::string readGroup;
	bool isSecondMate;
	TRecalibrationEMModelCovariateDefinition covariates;

	TRecalibrationEMModelDefinition(const std::string ReadGroup, const bool IsSecondMate){
		readGroup = ReadGroup;
		isSecondMate = IsSecondMate;
	};

	bool parseModel(const std::string modelString, std::string & error){
		return(covariates.parse(modelString, error));
	};
};

//--------------------------------------------------------------------
// TRecalibrationEMModel
//--------------------------------------------------------------------
class TRecalibrationEMModel{
private:
	TLog* logfile;

	//covariates
	TRecalibrationEMCovariateFunction_intercept intercept;
	std::vector<TRecalibrationEMCovariate*> covariates;
	std::vector<TRecalibrationEMCovariateFunction*> pointerToCovariateFunctions;

	//Newton Raphson Parameters
	//TODO: maybe split into class that can and cannot estimate?
	size_t _numParameters;
	double _Q, _oldQ;
	arma::mat Jacobian;
	arma::vec F;
	arma::mat JxF;
	TRecalibrationEMFirstDerivatives firstDerivatives;
	TRecalibrationEMSecondDerivatives secondDerivatives;
	unsigned int _numSitesAdded;
	bool _NRconverged;
	bool _NRStepAccepted;

	void _createCovariates(TRecalibrationEMModelCovariateDefinition & covariateMap, TRecalibrationEMDataTable* dataTable);
	void _createCovariates(TRecalibrationEMModelCovariateDefinition & covariateMap);
	void _initializeDerivatives();
	double _calcEpsilon(const double eta);
	inline double _calcQ(const int & genotype, TRecalibrationEMReadData & data){
		double eps = getErrorRate(data);
		double B = 1.33333333333333333333 * data.D[genotype] - 1.0;
		double P_d_given_g_beta = B * eps - data.D[genotype] + 1.0;
		return log(P_d_given_g_beta);
	};

public:
	TRecalibrationEMModel(TRecalibrationEMModelCovariateDefinition & covariateMap, TLog* Logfile);
	TRecalibrationEMModel(TRecalibrationEMModelCovariateDefinition & covariateMap, TRecalibrationEMDataTable* dataTable, TLog* Logfile);
	~TRecalibrationEMModel();

	bool checkParameterRange(TRecalibrationEMDataTable* dataTable, std::string & error);
	bool checkParameterRange(std::vector<uint16_t> & Qualities, uint16_t maxPos, std::string & error);
	int numParameters(){ return _numParameters; };

	void setEMParamsToZero();
	void setQToZero();
	void addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype);
	void addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta);
	double curQ(){ return _Q; };
	void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJFirstDerivatives, const double & weightJSecondDerivatives);
	bool solveJxF();
	void proposeNewParameters(double & lambda);
	bool acceptProposedParametersBasedOnQ();

	double getSteepestGradient();
	void printJacobianToStdOut();
	void printFToStdOut();
	void printJxFToStdOut();

	double getErrorRate(const TBase & base);
	double getErrorRate(const TRecalibrationEMReadData & data);

	TRecalibrationEMModelCovariateDefinition getCovariateDefinition();
};

//--------------------------------------------------------------------
// TRecalibrationEMModels
// Object containing a vector of recal models
//--------------------------------------------------------------------
class TRecalibrationEMModels{
private:
	TReadGroups* readGroups;
	TReadGroupMap* readGroupMap;
	std::vector<TRecalibrationEMModel> models;
	unsigned int totNumParameters;
	TRecalibrationEMReadGroupIndex readGroupIndex;
	TLog* logfile;


	void _readRecalFile(std::string filename, std::vector<TRecalibrationEMModelDefinition> & modelDefs);

	void _createModelsFromString(const std::string & s);
	void _createModelsFromFile(std::string filename);

	void _writeParameters(TOutputFile*out, const std::string & readGroupName, const int & readGroup, bool isSecondMate);

public:
	TRecalibrationEMModels(TReadGroups* ReadGroups, TReadGroupMap* readGroupMap, TLog* Logfile);
	~TRecalibrationEMModels();

	//general functions to add and remove models
	bool parseModelString(const std::string & modelString, std::map<std::string, std::string> covariateFunctions, std::string & error);
	void removeModel(int readGroupId, bool isSecondMate);

	//add model for estimation: dataTable provided
	void addModel(const int readGroupId, const bool isSecondMate, TRecalibrationEMModelCovariateDefinition & covariates, TRecalibrationEMDataTable* dataTable);
	void addModelsFromFile(std::string filename, TRecalibrationEMDataTables* dataTables);

	//add model for recalibration: no dataTable provided
	void addModel(const int readGroupId, const bool isSecondMate, const TRecalibrationEMModelCovariateDefinition & covariateFunctions);
	void createModels(std::string string);

	//inline TRecalibrationEMModel* operator[](int index){ return &models[index]; };
	int numModels(){ return models.size(); };
	bool modelExists(int readGroupId, bool isSecondMate){ return readGroupIndex.inUse(readGroupId, isSecondMate); };

	inline double calcEpsilon(const TRecalibrationEMReadData & data){
		return models[ readGroupIndex.index(data) ].getErrorRate(data);
	};
	inline double getErrorRate(TBase & base){
		return models[ readGroupIndex.index(base.readGroup, base.isSecondMate) ].getErrorRate(base);
	};

	bool hasReadGroupsWithoutModel();
	void addNoRecalModelIfMissing();
	void reportReadGroupsNotUsed();
	void reportReadGroupsConsideredSingleEnd();
	void warningForMissingReadGroups();

	//function to estimate
	void setEMParamsToZero();
	void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJFirstDerivatives, const double & weightJSecondDerivatives);
	void setQToZero();
	void addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta);
	void addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype);
	double curQ();
	bool solveJxF();
	void proposeNewParameters(double lambda);
	void scaleParameters();
	unsigned int acceptProposedParametersBasedOnQ();
	double getSteepestGradient();

	void writeRecalFile(TOutputFile* out);

};


#endif /* TRECALIBRATIONEMMODEL_H_ */
