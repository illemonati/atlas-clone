/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMMODEL_H_
#define TRECALIBRATIONEMMODEL_H_

#include "TFile.h"
#include "TSequencingErrorCovariate.h"
#include "auxiliaryTools.h"
#include "TGenotypeData.h"

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
// TSequencingErrorRho
//--------------------------------------------------------------------
class TSequencingErrorRho{
private:
	std::array< std::array<double, 4>, 4 > rho;

public:
	TSequencingErrorRho();
	void fillBaseLikelihoods(const Base base, const double epsilon, TBaseData & baseLikelihoods);
};

//--------------------------------------------------------------------
// TSequencingErrorModel
//--------------------------------------------------------------------
class TSequencingErrorModel{
private:
	TLog* logfile;

	//parameters: coavraites and rho
	TSequencingErrorCovariateList _covariates;
	TSequencingErrorRho rho;

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
	double _calcEpsilon(const double eta) const;
	inline double _calcQ(const double & eps, const Base & genotype, TRecalibrationEMReadData & data) const{
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

	void fillBaseLikelihoods(const TBase & base, TBaseData & baseLikelihoods);
	TSequencingErrorCovariateDefinition getCovariateDefinition();

	double getErrorRate(const TBase & base);
	double getErrorRate(const TRecalibrationEMReadData & data); //REMOVE!
};


}; //end namespace

#endif /* TRECALIBRATIONEMMODEL_H_ */
