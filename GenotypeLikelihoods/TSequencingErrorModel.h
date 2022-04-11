/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMMODEL_H_
#define TRECALIBRATIONEMMODEL_H_

#include <stdint.h>
#include <array>
#include <memory>
#include <string>
#include <vector>
#include <armadillo>

#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TGenotypeData.h"
#include "TSequencingErrorCovariate.h"
#include "TSequencingErrorCovariateFunction.h"
#include "auxiliaryTools.h"
#include "probability.h"

namespace BAM { class TSequencedBase; }
namespace GenotypeLikelihoods { namespace RecalEstimatorTools { class TRecalDataTable; } }

namespace GenotypeLikelihoods {
namespace SequencingError {


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

	double operator()(genometools::Base from, const genometools::Base to) const noexcept { return rho[genometools::index(from)][genometools::index(to)]; }
	std::string getDefinition() const noexcept;

	// functions used to estimate
	void prepareEstimationFromEMWeights() noexcept { rho.fill({0., 0., 0., 0.}); }
	void addBaseForEstimation(genometools::Base base, const TBaseLikelihoods &EMWeights) noexcept;
	void estimate() noexcept;
};

//--------------------------------------------------------------------
// TCovariateDef
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
struct TCovariateDef {
	std::string covariate;
	std::string function;
};

//--------------------------------------------------------------------
// TModelDefinition
// class to store model definition. Used when parsing files
//--------------------------------------------------------------------
class TModelDefinition {
public:
	std::vector<TCovariateDef> covariates; 
	std::string intercept;
	TRho rho;

	TModelDefinition() = default;
	TModelDefinition(const std::string &covariateString, const std::string &rhoString);
};

//--------------------------------------------------------------------
// TModel
// pure abstract base class
//--------------------------------------------------------------------
class TModel {
public:
	virtual ~TModel()                                                                                    = default;
	virtual bool estimatable() const noexcept                                                            = 0;
	virtual bool recalibrates() const noexcept                                                           = 0;
	virtual coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const noexcept          = 0;
	virtual genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const noexcept = 0;
	virtual void fillBaseLikelihoods(const BAM::TSequencedBase &base,
					 TBaseLikelihoods &baseLikelihoods) const noexcept                   = 0;
	virtual void simulate(BAM::TSequencedBase &base) const noexcept = 0;
	virtual std::string getCovariateDefinition() const noexcept                                                   = 0;
	virtual std::string getRhoDefinition() const noexcept                                                         = 0;
};

//------------------------------------------------
// TModelNoRecal
//------------------------------------------------
class TModelNoRecal : public TModel {
public:
	virtual bool estimatable() const noexcept override { return false; };
	virtual bool recalibrates() const noexcept override { return false; };

	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const noexcept override;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const noexcept override;
	void fillBaseLikelihoods(const BAM::TSequencedBase &base, TBaseLikelihoods &baseLikelihoods) const noexcept override;
	virtual void simulate(BAM::TSequencedBase &base) const noexcept override;

	virtual std::string getCovariateDefinition() const noexcept override { return "-"; };
	virtual std::string getRhoDefinition() const noexcept override { return "-"; };
};

//------------------------------------------------
// TModelRecal
//------------------------------------------------
class TModelRecal : public TModel {
private:
	struct TCovariateModel {
		std::unique_ptr<TCovariate> covariate;
		std::unique_ptr<TFunction> function;
		TCovariateModel(TCovariate* cov, TFunction* fn) : covariate(cov), function(fn) {}
	};
	TRho _rho;
	TIntercept _intercept;
	std::vector<TCovariateModel> _covariates;
	std::vector<TFunction *> _functions; // non-owning
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
	coretools::Probability _calcErrorRate(const BAM::TSequencedBase &base) const noexcept;

public:
	TModelRecal(const TModelDefinition &modelDef);
	TModelRecal(const TModelDefinition &modelDef, const RecalEstimatorTools::TRecalDataTable &DataTable);

	bool estimatable() const noexcept override { return true; };
	bool recalibrates() const noexcept override { return true; };
	std::string getCovariateDefinition() const noexcept override;
	std::string getRhoDefinition() const noexcept override { return _rho.getDefinition(); };
	TModelDefinition getModelDefinition() const;

	// get error rates
	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const noexcept override;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const noexcept override;
	void fillBaseLikelihoods(const BAM::TSequencedBase &base, TBaseLikelihoods &baseLikelihoods) const noexcept override;
	virtual void simulate(BAM::TSequencedBase &base) const noexcept override;

	// functions to estimate
	std::string checkParameterRange(RecalEstimatorTools::TRecalDataTable &DataTable) const;
	uint16_t numParameters() const noexcept { return _numParameters; }

	// functions to estimate rho
	void prepareRhoEstimationFromEMWeights() noexcept;
	void addBaseForRhoEstimation(BAM::TSequencedBase &base, const TBaseLikelihoods &EMWeights) noexcept;
	void estimateRho() noexcept;

	// functions to estimate betas
	void setNewtonRaphsonParamsToZero();
	void setQToZero() noexcept;
	void addToQ(const BAM::TSequencedBase &base, const TBaseLikelihoods &EM_weights_bbar_given_d);
	double curQ() const noexcept { return _Q; }
	void addToFandJacobian(const BAM::TSequencedBase &base, const TBaseLikelihoods &EM_weights_bbar_given_d);
	bool solveJxF();
	void proposeNewParameters(double &lambda);
	bool acceptProposedParametersBasedOnQ();
	void adjustParametersPostEstimation();
	double getSteepestGradient() const noexcept;

	const auto & Jacobian() const noexcept {return _Jacobian;}
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMMODEL_H_ */
