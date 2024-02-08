/*
 * TEpsilon.cpp
 *
 *  Created on: Jul 19, 2022
 *      Author: Andreas
 */

#include "TEpsilon.h"
#include "SequencingError/TCovariate.h"
#include "SequencingError/TFunction.h"
#include "SequencingError/TFunctions.h"
#include "coretools/Main/TError.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/splitters.h"
#include "coretools/Strings/toString.h"
#include <algorithm>
#include <iterator>
#include <memory>

namespace GenotypeLikelihoods {
namespace SequencingError {

using namespace coretools::str;

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

void TEpsilon::init(const RecalEstimatorTools::TRecalDataTable &DataTable) {
	_functions.init(DataTable);
	const size_t numParameters = _functions.numParameters();

	// prepare Newton-Raphson variables
	_Jacobian.resize(numParameters, numParameters);
	_F.resize(numParameters);
}

coretools::Probability TEpsilon::calcErrorRate(const BAM::TSequencedBase &base) const noexcept {
	return _functions.getEpsilon(base);
}

coretools::Probability TEpsilon::_calcErrorRate(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
											   std::vector<T2ndDerivative> &der2) const noexcept {
	// eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function
	return _functions.getEpsilon(base, der1, der2);
}


void TEpsilon::solveJxF() {
	// scale maxF #sites
	_maxF = std::max(_F.max(), -_F.min()) / _numSitesAdded;
	if (!solve(_JxF, _Jacobian, _F))
		UERROR("Issue solving JxF! This may be due to a lack of data. Consider adding more sites. Jacobian: ",
		       _Jacobian);

	// automatically reset
	_Jacobian.zeros();
	_F.zeros();
	_numSitesAdded = 0;
	_oldQ          = _Q;
}

void TEpsilon::propose(double lambda) {
	if (!_converged) {
		_functions.propose(lambda, _JxF);
		_Q = 0; // reset to recalculate
	}
}

bool TEpsilon::acceptOrReject() {
	_converged = _Q > _oldQ;
	if (!_converged) _functions.reject();
	return _converged;
}

void TEpsilon::adjust() {
	_Q         = 0.;
	_converged = false;
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
