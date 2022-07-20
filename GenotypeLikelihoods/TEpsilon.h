/*
 * TEpsilon.h
 *
 *  Created on: Jul 19, 2022
 *      Author: Andreas
 */

#ifndef TEPSILON_H_
#define TEPSILON_H_

#include <armadillo>
#include <memory>

#include "RecalEstimatorTools.h"
#include "TSequencingErrorCovariate.h"
#include "TSequencingErrorCovariateFunction.h"
#include "probability.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

struct TCovariateModel {
	std::unique_ptr<TCovariate> covariate;
	std::unique_ptr<TFunction> function;
	TCovariateModel(TCovariate *cov, TFunction *fn) : covariate(cov), function(fn) {}
};

class TEpsilon {
	TIntercept _intercept;
	std::vector<TCovariateModel> _covariates;
	size_t _numParameters;
	size_t _num1stDerivatives;
	size_t _num2ndDerivatives;

	double _Q    = 0.;
	double _oldQ = std::numeric_limits<double>::lowest();
	double _maxF = 0.;
	arma::mat _Jacobian;
	arma::vec _F;
	arma::mat _JxF;
	unsigned int _numSitesAdded = 0;

public:
	TEpsilon(const std::string& Def);

	void checkOrInit(const RecalEstimatorTools::TRecalDataTable &DataTable);

	coretools::Probability calcErrorRate(const BAM::TSequencedBase &base) const noexcept; 
	double Q() const noexcept {return _Q;};
	double maxF() const noexcept {return _maxF;};

	void addToEpsilon(const BAM::TSequencedBase &base, coretools::Probability P_g_I_d, coretools::Probability P_bbar_I_gd, bool updateJF=false);
	void solveJxF();
	void propose(double lambda);
	bool acceptOrReject();
	void adjust();

	std::string getDefinition() const noexcept;
};
} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMMODEL_H_ */
