/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMMODEL_H_
#define TRECALIBRATIONEMMODEL_H_

#include <armadillo>
#include <array>
#include <limits>
#include <memory>
#include <stdint.h>
#include <string>
#include <vector>

#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TGenotypeData.h"
#include "TSequencingErrorCovariate.h"
#include "TSequencingErrorCovariateFunction.h"
#include "TStrongArray.h"
#include "probability.h"
#include "devtools.h"

namespace BAM {
class TSequencedBase;
}
namespace GenotypeLikelihoods {
namespace RecalEstimatorTools {
class TRecalDataTable;
}
} // namespace GenotypeLikelihoods

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------------
// TRho
//--------------------------------------------------------------------
class TRho {
private:
	coretools::TStrongArray<TBaseProbabilities, genometools::Base> _rho{
		{coretools::TStrongArray<double, genometools::Base>{{0., 1. / 3, 1. / 3, 1. / 3}},
		 coretools::TStrongArray<double, genometools::Base>{{1. / 3, 0., 1. / 3, 1. / 3}},
		 coretools::TStrongArray<double, genometools::Base>{{1. / 3, 1. / 3, 0., 1. / 3}},
		 coretools::TStrongArray<double, genometools::Base>{{1. / 3, 1. / 3, 1. / 3, 0.}}}}; //[from][to]

	coretools::TStrongArray<coretools::TStrongArray<double, genometools::Base>, genometools::Base> _rhoSum{};
public:
	TRho() = default;
	TRho(const std::string &def);
	TRho(const TRho &other) = default;
	TRho &operator=(const TRho &other) = default;


	const TBaseProbabilities& operator[](genometools::Base from) const noexcept {
		return _rho[from];
	}
	std::string getDefinition() const noexcept;

	// functions used to estimate
	void add(genometools::Base base, coretools::Probability P_g_I_d, const TBaseProbabilities &P_bbar_I_d) noexcept;
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
	virtual TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept          = 0;
	virtual void simulate(BAM::TSequencedBase &base) const noexcept                                      = 0;
	virtual std::string getCovariateDefinition() const noexcept                                          = 0;
	virtual std::string getRhoDefinition() const noexcept                                                = 0;
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
	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept override;
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
		TCovariateModel(TCovariate *cov, TFunction *fn) : covariate(cov), function(fn) {}
	};
	TRho _rho;
	TIntercept _intercept;
	std::vector<TCovariateModel> _covariates;
	size_t _numParameters;
	size_t _num1stDerivatives;
	size_t _num2ndDerivatives;


	// Newton Raphson Parameters to estimate betas
	double _Q    = 0.;
	double _oldQ = std::numeric_limits<double>::lowest();
	double _maxF = 0.;
	arma::mat _Jacobian;
	arma::vec _F;
	arma::mat _JxF;
	unsigned int _numSitesAdded = 0;

	coretools::Probability _calcErrorRate(const BAM::TSequencedBase &base) const noexcept;

public:
	TModelRecal(const TModelDefinition &modelDef);

	bool estimatable() const noexcept override { return true; };
	bool recalibrates() const noexcept override { return true; };
	std::string getCovariateDefinition() const noexcept override;
	std::string getRhoDefinition() const noexcept override { return _rho.getDefinition(); };
	TModelDefinition getModelDefinition() const;

	// get error rates
	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const noexcept override;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const noexcept override;
	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept override;
	virtual void simulate(BAM::TSequencedBase &base) const noexcept override;

	// functions to estimate
	void checkOrInit(const RecalEstimatorTools::TRecalDataTable &DataTable) const;

	// functions to estimate rho
	void addToRho(const BAM::TSequencedBase &data, coretools::Probability P_g_I_d, const TBaseProbabilities &P_bbar_I_d) noexcept; 
	void estimateRho() noexcept;

	// functions to estimate betas
	void resetQ() noexcept {_oldQ = _Q; _Q = 0;};
	void addToQFJ(const BAM::TSequencedBase &base, coretools::Probability P_g_I_d, coretools::Probability P_bbar_I_gd, bool update=false);
	double curQ() const noexcept { return _Q; }
	void solveJxF();
	void proposeNewParameters(double lambda);
	bool acceptProposedParametersBasedOnQ();
	void adjustParametersPostEstimation();
	double maxF() const noexcept {return _maxF;};
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMMODEL_H_ */
