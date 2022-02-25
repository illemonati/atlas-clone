/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMMODEL_H_
#define TRECALIBRATIONEMMODEL_H_

#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TFile.h"
#include "TGenotypeData.h"
#include "TSequencingErrorCovariate.h"
#include "auxiliaryTools.h"
#include <memory>
#include <string>
#include <vector>

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------------
// TCovariateDef
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
struct TCovariateDef {
	std::string covariate;
	std::string function;
};

class TCovariateDefinition {
private:
	std::vector<TCovariateDef> _covariateFunctions; //<covariate, function>
	std::string _intercept = "";

public:
	TCovariateDefinition() = default;
	TCovariateDefinition(const std::string &modelString);

	void setIntercept(const double Intercept);
	void addCovariate(const std::string covariate, const std::string function);
	size_t size() const { return _covariateFunctions.size(); };
	const std::string &intercept() const { return _intercept; };
	auto begin() noexcept { return _covariateFunctions.begin(); };
	auto end() noexcept { return _covariateFunctions.end(); };
	auto cbegin() const noexcept { return _covariateFunctions.cbegin(); };
	auto cend() const noexcept { return _covariateFunctions.cend(); };
	std::string getModelString() const;
};

//--------------------------------------------------------------------
// TRho
//--------------------------------------------------------------------
class TRho {
private:
	std::array<std::array<double, 4>, 4> rho = {{{0., 1. / 3, 1. / 3, 1. / 3},
						     {1. / 3, 0., 1. / 3, 1. / 3},
						     {1. / 3, 1. / 3, 0., 1. / 3},
	                                             {1. / 3, 1. / 3, 1. / 3, 0.}}}; //[from][to]
public:
	TRho() = default;
	TRho(const std::string &def);
	TRho(const TRho &other) = default;
	TRho &operator=(const TRho &other) = default;

	double operator()(const uint8_t &from, const uint8_t &to) const noexcept { return rho[from][to]; }
public:
	std::string getDefinition() const noexcept;
	void fillBaseLikelihoods(const genometools::Base base, const coretools::Probability &epsilon,
				 TBaseLikelihoods &baseLikelihoods) const noexcept;

	// functions used to estimate
	void prepareEstimationFromEMWeights() noexcept { rho.fill({0., 0., 0., 0.}); }
	void addBaseForEstimation(const genometools::Base &base, const TBaseLikelihoods &EMWeights) noexcept;
	void estimate() noexcept;
};

//--------------------------------------------------------------------
// TModelDefinition
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
class TModelDefinition {
public:
	TCovariateDefinition covariates;
	TRho rho;

	TModelDefinition() = default;
	TModelDefinition(const std::string &covariateString, const std::string &rhoString)
		: covariates(covariateString), rho(rhoString) {}
};

//--------------------------------------------------------------------
// TModel
// pure abstract base class
//--------------------------------------------------------------------
class TModel {
public:
	virtual ~TModel()                                                                           = default;
	virtual bool estimatable() const noexcept                                                   = 0;
	virtual bool recalibrates() const noexcept                                                  = 0;
	virtual coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const          = 0;
	virtual genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const = 0;
	virtual void fillBaseLikelihoods(const BAM::TSequencedBase &base, TBaseLikelihoods &baseLikelihoods) const = 0;
	virtual void simulate(BAM::TSequencedBase &base, coretools::TRandomGenerator &RandomGenerator) const       = 0;
	virtual std::string getCovariateDefinition() const noexcept                                                = 0;
	virtual std::string getRhoDefinition() const noexcept                                                      = 0;
};

//------------------------------------------------
// TModelNoRecal
//------------------------------------------------
class TModelNoRecal : public TModel {
public:
	virtual bool estimatable() const noexcept override { return false; };
	virtual bool recalibrates() const noexcept override { return false; };

	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const override;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const override;
	void fillBaseLikelihoods(const BAM::TSequencedBase &base, TBaseLikelihoods &baseLikelihoods) const override;
	virtual void simulate(BAM::TSequencedBase &base, coretools::TRandomGenerator &RandomGenerator) const override;

	virtual std::string getCovariateDefinition() const noexcept override { return "-"; };
	virtual std::string getRhoDefinition() const noexcept override { return "-"; };
};

//------------------------------------------------
// TModelRecal
//------------------------------------------------
class TModelRecal : public TModel {
private:
	TRho _rho;
	TCovariateFunction_intercept _intercept;
	std::vector<std::unique_ptr<TCovariate>> _covariates;
	std::vector<TCovariateFunction *> _functions; // non-owning
	uint16_t _numParameters;

	// Newton Raphson Parameters to estimate betas
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
	coretools::Probability _calcErrorRate(const BAM::TSequencedBase &base) const;

public:
	TModelRecal(const TModelDefinition &modelDef);
	TModelRecal(const TModelDefinition &modelDef, const RecalEstimatorTools::TRecalDataTable &DataTable);

	bool estimatable() const noexcept override { return true; };
	bool recalibrates() const noexcept override { return true; };
	std::string getCovariateDefinition() const noexcept override;
	std::string getRhoDefinition() const noexcept override { return _rho.getDefinition(); };
	TModelDefinition getModelDefinition() const;

	// get error rates
	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const override;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const override;
	void fillBaseLikelihoods(const BAM::TSequencedBase &base, TBaseLikelihoods &baseLikelihoods) const override;
	virtual void simulate(BAM::TSequencedBase &base, coretools::TRandomGenerator &RandomGenerator) const override{};

	// functions to estimate
	bool checkParameterRange(RecalEstimatorTools::TRecalDataTable &DataTable, std::string &error);
	uint16_t numParameters() { return _numParameters; };

	// functions to estimate rho
	void prepareRhoEstimationFromEMWeights();
	void addBaseForRhoEstimation(BAM::TSequencedBase &base, const TBaseLikelihoods &EMWeights);
	void estimateRho();

	// functions to estimate betas
	void setNewtonRaphsonParamsToZero();
	void setQToZero();
	void addToQ(const BAM::TSequencedBase &base, const TBaseLikelihoods &EM_weights_bbar_given_d);
	double curQ() { return _Q; };
	void addToFandJacobian(const BAM::TSequencedBase &base, const TBaseLikelihoods &EM_weights_bbar_given_d);
	bool solveJxF();
	void proposeNewParameters(double &lambda);
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
