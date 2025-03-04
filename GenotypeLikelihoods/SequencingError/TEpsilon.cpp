/*
 * TEpsilon.cpp
 *
 *  Created on: Jul 19, 2022
 *      Author: Andreas
 */

#include "TEpsilon.h"
#include "armadillo"

namespace GenotypeLikelihoods {
namespace SequencingError {

TEpsilon::TEpsilon(std::string_view Def) : _functions(Def) {
	const size_t numParameters = _functions.numParameters();

	// prepare Newton-Raphson variables
	_Jacobian.resize(numParameters, numParameters);
	_F.resize(numParameters);
}

TEpsilon::TEpsilon(const BAM::RGInfo::TInfo &Info) : _functions(Info) {
	const size_t numParameters = _functions.numParameters();

	// prepare Newton-Raphson variables
	_Jacobian.resize(numParameters, numParameters);
	_F.resize(numParameters);
}

void TEpsilon::init(const RecalEstimatorTools::TRecalDataTable &DataTable, size_t MinData) {
	_functions.init(DataTable, MinData);
	const size_t numParameters = _functions.numParameters();

	// prepare Newton-Raphson variables
	_Jacobian.resize(numParameters, numParameters);
	_F.resize(numParameters);
}

coretools::Probability TEpsilon::calcErrorRate(const BAM::TSequencedData &data) const noexcept {
	return _functions.getEpsilon(data);
}

coretools::Probability TEpsilon::_calcErrorRate(const BAM::TSequencedData &data, std::vector<T1stDerivative> &der1,
											   std::vector<T2ndDerivative> &der2) const noexcept {
	// eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function
	return _functions.getEpsilon(data, der1, der2);
}


void TEpsilon::solveJxF() {
	// scale maxF #sites
	_Jacobian = arma::symmatu(_Jacobian);
	_maxF = std::max(std::abs(_F.max()), std::abs(_F.min())) / _numSitesAdded;
	if (!solve(_JxF, _Jacobian, _F))
		UERROR("Issue solving JxF! This may be due to a lack of data. Consider adding more sites. Jacobian: ",
		       _Jacobian);
	_maxJxF = std::max(std::abs(_JxF.max()), std::abs(_JxF.min()));

	// automatically reset
	_Jacobian.zeros();
	_F.zeros();
	_numSitesAdded = 0;
	_oldQ          = _Q;
}

void TEpsilon::propose(double lambda) {
	if (!_accepted) {
		_functions.propose(lambda, _JxF);
		_Q = 0; // reset to recalculate
	}
}

bool TEpsilon::acceptOrReject() {
	_accepted = _Q > _oldQ;
	if (!_accepted) _functions.reject();
	return _accepted;
}

void TEpsilon::adjust() {
	_Q         = 0.;
	_accepted = false;
	_functions.adjust();
}

void TEpsilon::log() const {
	_functions.log();
}

BAM::RGInfo::TInfo TEpsilon::info() const{
	return _functions.info();
}

} // namespace SequencingError
} // namespace GenotypeLikelihoods
