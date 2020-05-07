/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMMODEL_H_
#define TRECALIBRATIONEMMODEL_H_

#include "TFile.h"
#include "../GenotypeLikelihoods/TSequencingErrorCovariate.h"
#include "auxiliaryTools.h"

namespace GenotypeLikelihoods{

//--------------------------------------------------------------------
// TSequencingErrorCovariateDef
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
struct TSequencingErrorCovariateDef{
	std::string covariate;
	std::string function;

	TSequencingErrorCovariateDef(const std::string Covariate, const std::string Function){
		covariate = Covariate;
		function = Function;
	};
};

typedef std::vector<TSequencingErrorCovariateDef>::iterator TRecalibrationEMModelCovariateDefinitionIterator;
class TSequencingErrorCovariateDefinition{
private:
	std::vector<TSequencingErrorCovariateDef> covariateFunctions;  //<covariate, function>

public:
	std::string intercept;

	TSequencingErrorCovariateDefinition(){
		intercept = "";
	};
	TSequencingErrorCovariateDefinition(const std::string modelString, std::string & error){
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
// TSequencingErrorModelDefinition
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
class TSequencingErrorModelDefinition{
public:
	uint16_t readGroupId;
	bool isSecondMate;
	TSequencingErrorCovariateDefinition covariates;

	TSequencingErrorModelDefinition(const uint16_t ReadGroupId, const bool IsSecondMate){
		readGroupId = ReadGroupId;
		isSecondMate = IsSecondMate;
	};

	bool parseModel(const std::string modelString, std::string & error){
		return(covariates.parse(modelString, error));
	};
};

//--------------------------------------------------------------------
// TSequencingErrorCovariateList
//--------------------------------------------------------------------
class TSequencingErrorCovariateList{
private:
	void _storePointersToCovariateFunctions();
	void _clear();

public:
	uint16_t numParameters;
	TSequencingErrorCovariateFunction_intercept intercept;
	std::vector< TSequencingErrorCovariate* > covariates;
	std::vector< TSequencingErrorCovariateFunction* > pointerToCovariateFunctions;

	TSequencingErrorCovariateList();
	~TSequencingErrorCovariateList();
	TSequencingErrorCovariateList(TSequencingErrorCovariateList&& other);
	TSequencingErrorCovariateList& operator=(TSequencingErrorCovariateList&& other);

	void _createCovariatesAndIntercept(TSequencingErrorCovariateDefinition & covariateMap, TRecalibrationEMDataTable* dataTable);
	void _createCovariatesAndIntercept(TSequencingErrorCovariateDefinition & covariateMap);
	TSequencingErrorCovariateDefinition getCovariateDefinition();
};

//--------------------------------------------------------------------
// TSequencingErrorModel
//--------------------------------------------------------------------
class TSequencingErrorModel{
private:
	TLog* logfile;

	//covariates
	TSequencingErrorCovariateList _covariates;

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
	TSequencingErrorModel(TSequencingErrorCovariateDefinition & covariateMap, TLog* Logfile);
	TSequencingErrorModel(TSequencingErrorCovariateDefinition & covariateMap, TRecalibrationEMDataTable* dataTable, TLog* Logfile);

	bool checkParameterRange(TRecalibrationEMDataTable* dataTable, std::string & error);
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

	double getErrorRate(const TBaseData & base);
	double getErrorRate(const TRecalibrationEMReadData & data);

	TSequencingErrorCovariateDefinition getCovariateDefinition();
};

//--------------------------------------------------------------------
// TSequencingErrorModels
// Object containing a vector of recal models
//--------------------------------------------------------------------
class TSequencingErrorModels{
private:
	TLog* logfile;
	TReadGroups* readGroups;
	TReadGroupMap* readGroupMap;
	std::vector<TSequencingErrorModel> models;
	unsigned int totNumParameters;
	TRecalibrationEMReadGroupIndex readGroupIndex;

	void _readRecalFile(const std::string filename, std::vector<TSequencingErrorModelDefinition> & modelDefs);

	void _createModelsFromString(const std::string & s);
	void _createModelsFromFile(std::string filename);

	void _writeParameters(TOutputFile & out, const std::string & readGroupName, const int & readGroup, bool isSecondMate);

public:
	TSequencingErrorModels(TReadGroups* ReadGroups, TReadGroupMap* readGroupMap, TLog* Logfile);

	//general functions to add and remove models
	//bool parseModelString(const std::string & modelString, std::map<std::string, std::string> covariateFunctions, std::string & error);
	void removeModel(int readGroupId, bool isSecondMate);

	//add model for estimation: dataTable provided
	void addModel(const uint16_t readGroupId, const bool isSecondMate, TSequencingErrorCovariateDefinition & covariates, TRecalibrationEMDataTable* dataTable);
	void addModelsFromFile(std::string filename, TRecalibrationEMDataTables* dataTables);

	//add model for recalibration: no dataTable provided
	void addModel(const uint16_t readGroupId, const bool isSecondMate, TSequencingErrorCovariateDefinition & covariateFunctions);
	void createModels(std::string string);

	//inline TRecalibrationEMModel* operator[](int index){ return &models[index]; };
	int numModels(){ return models.size(); };
	bool modelExists(uint16_t readGroupId, bool isSecondMate){ return readGroupIndex.inUse(readGroupId, isSecondMate); };
	bool modelExists(TSequencingErrorModelDefinition & def){ return readGroupIndex.inUse(def.readGroupId, def.isSecondMate); };

	inline double calcEpsilon(const TRecalibrationEMReadData & data){
		return models[ readGroupIndex.index(data) ].getErrorRate(data);
	};
	inline double getErrorRate(TBaseData & base){
		return models[ readGroupIndex.index(base.readGroup, base.isSecondMate()) ].getErrorRate(base);
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

}; //end namespace

#endif /* TRECALIBRATIONEMMODEL_H_ */
