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
struct TRecalibrationEMModelCovariateDef{
	std::string covariate;
	std::string function;

	TRecalibrationEMModelCovariateDef(const std::string Covariate, const std::string Function){
		covariate = Covariate;
		function = Function;
	};
};

typedef std::vector<TRecalibrationEMModelCovariateDef>::iterator TRecalibrationEMModelCovariateDefinitionIterator;
class TRecalibrationEMModelCovariateDefinition{
private:
	std::vector<TRecalibrationEMModelCovariateDef> covariateFunctions;  //<covariate, function>

public:
	std::string intercept;

	TRecalibrationEMModelCovariateDefinition(){
		intercept = "";
	};
	TRecalibrationEMModelCovariateDefinition(const std::string modelString, std::string & error){
		parse(modelString, error);
	};

	bool parse(const std::string & modelString, std::string & error);
	void setIntercept(const double Intercept);
	void addCovariate(const std::string covariate, const std::string function);
	size_t size(){ return covariateFunctions.size(); };
	TRecalibrationEMModelCovariateDefinitionIterator begin(){ return covariateFunctions.begin(); };
	TRecalibrationEMModelCovariateDefinitionIterator end(){ return covariateFunctions.end(); };
	std::string getModelString();
};

//--------------------------------------------------------------------
// TRecalibrationEMModelDefinition
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
class TRecalibrationEMModelDefinition{
public:
	uint16_t readGroupId;
	bool isSecondMate;
	TRecalibrationEMModelCovariateDefinition covariates;

	TRecalibrationEMModelDefinition(const uint16_t ReadGroupId, const bool IsSecondMate){
		readGroupId = ReadGroupId;
		isSecondMate = IsSecondMate;
	};

	bool parseModel(const std::string modelString, std::string & error){
		return(covariates.parse(modelString, error));
	};
};

//--------------------------------------------------------------------
// TRecalibrationEMModelCovariateList
//--------------------------------------------------------------------
class TRecalibrationEMModelCovariateList{
private:
	void _storePointersToCovariateFunctions();
	void _clear();

public:
	uint16_t numParameters;
	TRecalibrationEMCovariateFunction_intercept intercept;
	std::vector< TRecalibrationEMCovariate* > covariates;
	std::vector< TRecalibrationEMCovariateFunction* > pointerToCovariateFunctions;

	TRecalibrationEMModelCovariateList();
	~TRecalibrationEMModelCovariateList();
	TRecalibrationEMModelCovariateList(TRecalibrationEMModelCovariateList&& other);
	TRecalibrationEMModelCovariateList& operator=(TRecalibrationEMModelCovariateList&& other);

	void _createCovariatesAndIntercept(TRecalibrationEMModelCovariateDefinition & covariateMap, TRecalibrationEMDataTable* dataTable);
	void _createCovariatesAndIntercept(TRecalibrationEMModelCovariateDefinition & covariateMap);
	TRecalibrationEMModelCovariateDefinition getCovariateDefinition();
};

//--------------------------------------------------------------------
// TRecalibrationEMModel
//--------------------------------------------------------------------
class TRecalibrationEMModel{
private:
	TLog* logfile;

	//covariates
	TRecalibrationEMModelCovariateList _covariates;

	//Newton Raphson Parameters
	double _Q, _oldQ;
	arma::mat _Jacobian;
	arma::vec _F;
	arma::mat _JxF;
	TRecalibrationEMFirstDerivatives _firstDerivatives;
	TRecalibrationEMSecondDerivatives _secondDerivatives;
	unsigned int _numSitesAdded;
	bool _NRconverged;
	bool _NRStepAccepted;

	void _initializeDerivatives();
	double _calcEpsilon(const double eta);
	inline double _calcQ(const double & eps, const Base & genotype, TRecalibrationEMReadData & data){
		double B = 1.33333333333333333333 * data.D[genotype] - 1.0;
		double P_d_given_g_beta = B * eps - data.D[genotype] + 1.0;
		return log(P_d_given_g_beta);
	};

public:
	TRecalibrationEMModel(TRecalibrationEMModelCovariateDefinition & covariateMap, TLog* Logfile);
	TRecalibrationEMModel(TRecalibrationEMModelCovariateDefinition & covariateMap, TRecalibrationEMDataTable* dataTable, TLog* Logfile);

	bool checkParameterRange(TRecalibrationEMDataTable* dataTable, std::string & error);
	bool checkParameterRange(std::vector<uint16_t> & Qualities, uint16_t maxPos, std::string & error);
	uint16_t numParameters(){ return _covariates.numParameters; };

	void setEMParamsToZero();
	void setQToZero();
	void addToQ(TRecalibrationEMReadData & data, const Base & knownGenotype);
	void addToQ(TRecalibrationEMReadData & data, double* P_g_given_d_oldBeta);
	double curQ(){ return _Q; };
	void addToFandJacobian(const TRecalibrationEMReadData & data, const double & weightF, const double & weightJacobian);
	bool solveJxF();
	void proposeNewParameters(double & lambda);
	bool acceptProposedParametersBasedOnQ();
	void adjustParametersPostEstimation();
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
	TLog* logfile;
	TReadGroups* readGroups;
	TReadGroupMap* readGroupMap;
	std::vector<TRecalibrationEMModel> models;
	unsigned int totNumParameters;
	TRecalibrationEMReadGroupIndex readGroupIndex;

	void _readRecalFile(const std::string filename, std::vector<TRecalibrationEMModelDefinition> & modelDefs);

	void _createModelsFromString(const std::string & s);
	void _createModelsFromFile(std::string filename);

	void _writeParameters(TOutputFile & out, const std::string & readGroupName, const int & readGroup, bool isSecondMate);

public:
	TRecalibrationEMModels(TReadGroups* ReadGroups, TReadGroupMap* readGroupMap, TLog* Logfile);

	//general functions to add and remove models
	//bool parseModelString(const std::string & modelString, std::map<std::string, std::string> covariateFunctions, std::string & error);
	void removeModel(int readGroupId, bool isSecondMate);

	//add model for estimation: dataTable provided
	void addModel(const uint16_t readGroupId, const bool isSecondMate, TRecalibrationEMModelCovariateDefinition & covariates, TRecalibrationEMDataTable* dataTable);
	void addModelsFromFile(std::string filename, TRecalibrationEMDataTables* dataTables);

	//add model for recalibration: no dataTable provided
	void addModel(const uint16_t readGroupId, const bool isSecondMate, TRecalibrationEMModelCovariateDefinition & covariateFunctions);
	void createModels(std::string string);

	//inline TRecalibrationEMModel* operator[](int index){ return &models[index]; };
	int numModels(){ return models.size(); };
	bool modelExists(uint16_t readGroupId, bool isSecondMate){ return readGroupIndex.inUse(readGroupId, isSecondMate); };
	bool modelExists(TRecalibrationEMModelDefinition & def){ return readGroupIndex.inUse(def.readGroupId, def.isSecondMate); };

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

	void writeRecalFile(const std::string filename);

};


#endif /* TRECALIBRATIONEMMODEL_H_ */
