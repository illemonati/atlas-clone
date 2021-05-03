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
	std::vector<TSequencingErrorCovariateDef> _covariateFunctions;  //<covariate, function>
	std::string _intercept;

public:
	TSequencingErrorCovariateDefinition(){
		clear();
	};
	TSequencingErrorCovariateDefinition(const std::string modelString, std::string & error){
		parse(modelString, error);
	};

	void clear();
	bool parse(const std::string & modelString, std::string & error);
	void setIntercept(const double Intercept);
	void addCovariate(const std::string covariate, const std::string function);
	size_t size(){ return _covariateFunctions.size(); };
	const std::string& intercept() const { return _intercept; };
	TRecalibrationEMModelCovariateDefinitionIterator begin(){ return _covariateFunctions.begin(); };
	TRecalibrationEMModelCovariateDefinitionIterator end(){ return _covariateFunctions.end(); };
	std::string getModelString();
};

//--------------------------------------------------------------------
// TSequencingErrorRhoStorage
//--------------------------------------------------------------------
class TSequencingErrorRhoStorage{
protected:
	std::array< std::array<double, 4>, 4 > rho; //[from][to]

public:
	TSequencingErrorRhoStorage();
	TSequencingErrorRhoStorage(const TSequencingErrorRhoStorage & other);
	TSequencingErrorRhoStorage(const std::string & def, std::string & error);
	void operator=(const TSequencingErrorRhoStorage & other);
	double operator()(const uint8_t & from, const uint8_t & to);
	void reset();
	bool set(const std::string & def, std::string & error);
	void set(const std::string & def);
	std::string getDefinition() const;
};

//--------------------------------------------------------------------
// TSequencingErrorRho
//--------------------------------------------------------------------
class TSequencingErrorRho:public TSequencingErrorRhoStorage{
public:
	void operator=(const TSequencingErrorRhoStorage & other);
	void fillBaseLikelihoods(const Base base, const double epsilon, TBaseData & baseLikelihoods) const;

	//functions used to estimate
	void prepareEstimationFromEMWeights();
	void addBaseForEstimation(const Base & base, const TBaseData & EMWeights);
	void estimate();
};

//--------------------------------------------------------------------
// TSequencingErrorModelDefinition
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
class TSequencingErrorModelDefinition{
public:
	TSequencingErrorCovariateDefinition covariates;
	TSequencingErrorRhoStorage rho;

	TSequencingErrorModelDefinition() = default;

	bool parseCovariates(const std::string & covariateString, std::string & error){
		return(covariates.parse(covariateString, error));
	};

	bool parseRho(const std::string & rhoString, std::string & error){
		return rho.set(rhoString, error);
	};

	bool parse(const std::string & covariateString, const std::string & rhoString, std::string & error){
		if(!parseCovariates(covariateString, error)){
			return false;
		}
		return parseRho(rhoString, error);
	};
};

//--------------------------------------------------------------------
// TSequencingErrorCovariateList
// class to store all covariates of an error model
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

	void createCovariatesAndIntercept(TSequencingErrorCovariateDefinition & covariateMap, TRecalibrationEMDataTable* dataTable);
	void createCovariatesAndIntercept(TSequencingErrorCovariateDefinition & covariateMap);
	TSequencingErrorCovariateDefinition getCovariateDefinition() const;
};

//--------------------------------------------------------------------
// TSequencingErrorModel
// pure abstract base class
//--------------------------------------------------------------------
class TSequencingErrorModel{
protected:
	TSequencingErrorRho _rho;

public:
	TSequencingErrorModel() = default;
	virtual ~TSequencingErrorModel(){};

	virtual bool estimatable() const { return false; };
	virtual bool recalibrates() const { return false; };

	virtual double getErrorRate(const BAM::TBase & base, const BAM::TQualityMap & qualMap) const = 0;
	virtual uint8_t getPhredInt(const BAM::TBase & base, const BAM::TQualityMap & qualMap) const = 0;
	virtual void fillBaseLikelihoods(const BAM::TBase & base, const BAM::TQualityMap & qualMap, TBaseData & baseLikelihoods) const = 0;

	virtual std::string getCovariateDefinition() const { return "-"; };
	virtual std::string getRhoDefinition() const { return "-"; };
};

//------------------------------------------------
// TSequencingErrorModelNoRecal
//------------------------------------------------
class TSequencingErrorModelNoRecal:public TSequencingErrorModel{
private:

public:
	TSequencingErrorModelNoRecal() = default;
	~TSequencingErrorModelNoRecal() = default;

	double getErrorRate(const BAM::TBase & base, const BAM::TQualityMap & qualMap) const override;
	uint8_t getPhredInt(const BAM::TBase & base, const BAM::TQualityMap & qualMap) const override;
	void fillBaseLikelihoods(const BAM::TBase & base, const BAM::TQualityMap & qualMap, TBaseData & baseLikelihoods) const override;
};

//------------------------------------------------
// TSequencingErrorModelRecal
//------------------------------------------------
class TSequencingErrorModelRecal:public TSequencingErrorModel{
private:
	TLog* _logfile;

	//parameters: covarites and rho
	TSequencingErrorCovariateList _covariates;

	//Newton Raphson Parameters to estimate betas
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

public:
	TSequencingErrorModelRecal(TSequencingErrorModelDefinition & modelDef, TLog* Logfile);
	TSequencingErrorModelRecal(TSequencingErrorModelDefinition & modelDef, TRecalibrationEMDataTable* dataTable, TLog* Logfile);

	bool estimatable() const override { return true; };
	bool recalibrates() const override { return true; };

	bool checkParameterRange(TRecalibrationEMDataTable* dataTable, std::string & error);
	uint16_t numParameters(){ return _covariates.numParameters; };

	//functions to estimate rho
	void prepareRhoEstimationFromEMWeights();
	void addBaseForRhoEstimation(BAM::TBase & base, const TBaseData & EMWeights);
	void estimateRho();

	//functions to estimate betas
	void setNewtonRaphsonParamsToZero();
	void setQToZero();
	void addToQ(const BAM::TBase & base, const TBaseData & EM_weights_bbar_given_d);
	double curQ(){ return _Q; };
	void addToFandJacobian(const BAM::TBase & base, const TBaseData & EM_weights_bbar_given_d);
	bool solveJxF();
	void proposeNewParameters(double & lambda);
	bool acceptProposedParametersBasedOnQ();
	void adjustParametersPostEstimation();
	double getSteepestGradient();
	void printJacobianToStdOut();
	void printFToStdOut();
	void printJxFToStdOut();

	void fillBaseLikelihoods(const BAM::TBase & base, TBaseData & baseLikelihoods) const;
	std::string getCovariateDefinition() const override;
	std::string getRhoDefinition() const override;

	double getErrorRate(const BAM::TBase & base) const;
};

}; //end namespace

#endif /* TRECALIBRATIONEMMODEL_H_ */
