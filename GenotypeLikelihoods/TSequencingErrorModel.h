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
#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include <vector>
#include <string>

namespace GenotypeLikelihoods{
namespace SequencingError {


//--------------------------------------------------------------------
// TCovariateDef
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
struct TCovariateDef{
	std::string covariate;
	std::string function;

	TCovariateDef(const std::string Covariate, const std::string Function){
		covariate = Covariate;
		function = Function;
	};
};

class TCovariateDefinition{
private:
	std::vector<TCovariateDef> _covariateFunctions;  //<covariate, function>
	std::string _intercept;

public:
	TCovariateDefinition(){
		clear();
	};
	TCovariateDefinition(const std::string modelString, std::string & error){
		parse(modelString, error);
	};

	void clear();
	bool parse(const std::string & modelString, std::string & error);
	void setIntercept(const double Intercept);
	void addCovariate(const std::string covariate, const std::string function);
	size_t size() const { return _covariateFunctions.size(); };
	const std::string& intercept() const { return _intercept; };
	std::vector<TCovariateDef>::iterator begin(){ return _covariateFunctions.begin(); };
	std::vector<TCovariateDef>::iterator end(){ return _covariateFunctions.end(); };
	std::vector<TCovariateDef>::const_iterator cbegin() const { return _covariateFunctions.cbegin(); };
	std::vector<TCovariateDef>::const_iterator cend() const { return _covariateFunctions.cend(); };
	std::string getModelString() const;
};

//--------------------------------------------------------------------
// TRhoStorage
//--------------------------------------------------------------------
class TRhoStorage{
protected:
	std::array< std::array<double, 4>, 4 > rho; //[from][to]

public:
	TRhoStorage();
	TRhoStorage(const TRhoStorage & other);
	TRhoStorage(const std::string & def, std::string & error);
	void operator=(const TRhoStorage & other);
	double operator()(const uint8_t & from, const uint8_t & to);
	void reset();
	bool set(const std::string & def, std::string & error);
	void set(const std::string & def);
	std::string getDefinition() const;
};

//--------------------------------------------------------------------
// TRho
//--------------------------------------------------------------------
class TRho:public TRhoStorage{
public:
	void operator=(const TRhoStorage & other);
	void fillBaseLikelihoods(const genometools::Base base, const coretools::Probability & epsilon, TBaseLikelihoods & baseLikelihoods) const;

	//functions used to estimate
	void prepareEstimationFromEMWeights();
	void addBaseForEstimation(const genometools::Base & base, const TBaseLikelihoods & EMWeights);
	void estimate();
};

//--------------------------------------------------------------------
// TModelDefinition
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
class TModelDefinition{
public:
	TCovariateDefinition covariates;
	TRhoStorage rho;

	TModelDefinition() = default;
	TModelDefinition(const std::string & covariateString, const std::string & rhoString, std::string & error){
		parse(covariateString, rhoString, error);
	};

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
// TCovariateList
// class to store all covariates of an error model
//--------------------------------------------------------------------
class TCovariateList{
private:
	void _storePointersToCovariateFunctions();
	void _clear();

public:
	uint16_t numParameters;
	TCovariateFunction_intercept intercept;
	std::vector< TCovariate* > covariates;
	std::vector< TCovariateFunction* > pointerToCovariateFunctions;

	TCovariateList();
	~TCovariateList();
	TCovariateList(TCovariateList&& other);
	TCovariateList& operator=(TCovariateList&& other);

	void createCovariatesAndIntercept(const TCovariateDefinition & covariateMap, const RecalEstimatorTools::TRecalDataTable & DataTable);
	void createCovariatesAndIntercept(const TCovariateDefinition & covariateMap);
	TCovariateDefinition getCovariateDefinition() const;
};

//--------------------------------------------------------------------
// TModel
// pure abstract base class
//--------------------------------------------------------------------
class TModel{
protected:
	TRho _rho;

public:
	TModel() = default;
	virtual ~TModel(){};

	virtual bool estimatable() const { return false; };
	virtual bool recalibrates() const { return false; };

	virtual coretools::Probability getErrorRate(const BAM::TSequencedBase & base) const = 0;
	virtual genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase & base) const = 0;
	virtual void fillBaseLikelihoods(const BAM::TSequencedBase & base, TBaseLikelihoods & baseLikelihoods) const = 0;

	virtual std::string getCovariateDefinition() const { return "-"; };
	virtual std::string getRhoDefinition() const { return "-"; };
};

//------------------------------------------------
// TModelNoRecal
//------------------------------------------------
class TModelNoRecal:public TModel{
private:

public:
	TModelNoRecal() = default;
	~TModelNoRecal() = default;

	coretools::Probability getErrorRate(const BAM::TSequencedBase & base) const override;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase & base) const override;
	void fillBaseLikelihoods(const BAM::TSequencedBase & base, TBaseLikelihoods & baseLikelihoods) const override;
};

//------------------------------------------------
// TModelRecal
//------------------------------------------------
class TModelRecal:public TModel{
private:
	//parameters: covarites and rho
	TCovariateList _covariates;

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
	coretools::Probability _calcEpsilon(double eta) const;
	coretools::Probability _calcErrorRate(const BAM::TSequencedBase & base) const;

public:
	TModelRecal(const TModelDefinition & modelDef);
	TModelRecal(const TModelDefinition & modelDef, const RecalEstimatorTools::TRecalDataTable & DataTable);

	bool estimatable() const override { return true; };
	bool recalibrates() const override { return true; };
	std::string getCovariateDefinition() const override { return _covariates.getCovariateDefinition().getModelString(); };
	std::string getRhoDefinition() const override { return _rho.getDefinition(); };
	TModelDefinition getModelDefinition() const;

	//get error rates
	coretools::Probability getErrorRate(const BAM::TSequencedBase & base) const override;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase & base) const override;
	void fillBaseLikelihoods(const BAM::TSequencedBase & base, TBaseLikelihoods & baseLikelihoods) const override;

	//functions to estimate
	bool checkParameterRange(RecalEstimatorTools::TRecalDataTable & DataTable, std::string & error);
	uint16_t numParameters(){ return _covariates.numParameters; };

	//functions to estimate rho
	void prepareRhoEstimationFromEMWeights();
	void addBaseForRhoEstimation(BAM::TSequencedBase & base, const TBaseLikelihoods & EMWeights);
	void estimateRho();

	//functions to estimate betas
	void setNewtonRaphsonParamsToZero();
	void setQToZero();
	void addToQ(const BAM::TSequencedBase & base, const TBaseLikelihoods & EM_weights_bbar_given_d);
	double curQ(){ return _Q; };
	void addToFandJacobian(const BAM::TSequencedBase & base, const TBaseLikelihoods & EM_weights_bbar_given_d);
	bool solveJxF();
	void proposeNewParameters(double & lambda);
	bool acceptProposedParametersBasedOnQ();
	void adjustParametersPostEstimation();
	double getSteepestGradient();
	void printJacobianToStdOut();
	void printFToStdOut();
	void printJxFToStdOut();
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMMODEL_H_ */
