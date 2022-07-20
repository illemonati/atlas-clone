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
#include "probability.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

class TFunction;

class TEpsilon {
	std::vector<std::unique_ptr<TFunction>> _functions;
	size_t _numParameters     = 0;
	size_t _num1stDerivatives = 0;
	size_t _num2ndDerivatives = 0;

	double _Q    = 0.;
	double _oldQ = std::numeric_limits<double>::lowest();
	double _maxF = 0.;
	arma::mat _Jacobian;
	arma::vec _F;
	arma::mat _JxF;
	unsigned int _numSitesAdded = 0;

public:
	TEpsilon(const std::string& Def);
	~TEpsilon();

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
