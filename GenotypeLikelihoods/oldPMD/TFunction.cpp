#include "TFunction.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"
#include <tuple>

namespace GenotypeLikelihoods::oldPMD {

using namespace coretools::str;
using coretools::Probability;
using coretools::instances::logfile;
using coretools::instances::parameters;

namespace impl {

template<typename T>
std::vector<T> parseParameters(std::string_view string) {
	// expect string of the form NAME[P1,P2,...]
	// extract P1, P2, ... as a vector of doubles
	std::vector<T> ps;
	if (string.find('[') == string.npos) return ps;

	auto tmp = readBefore(readAfter(string, '['), ']');
	if (tmp.empty()) return ps;

	TSplitter spl(readBefore(readAfter(string, '['), ']'), ',');

	for (auto s: spl) {
		ps.push_back(fromString<T, true>(strip(s)));
	}
	return ps;
}

std::array<double, 3> initialEstimatesOLS(const std::vector<Probability> &Empiric) {
	// fill vector y to fit using OLS
	const auto N = Empiric.size();
	arma::vec y(Empiric.size());
	double sumYSquared = 0.0;
	for (size_t p = 0; p < N; ++p) {
		y(p) = Empiric[p];
		sumYSquared += y(p) * y(p);
	}

	// some variables
	double gammaStep = 0.01;
	double gammaTmp  = -gammaStep + 0.00000001;
	double SSRold    = N;
	double SSRdiff   = -1.0;
	arma::mat betaHat;
	arma::mat X(N, 2);
	X.ones();

	// do until we get a small alpha
	while (fabs(gammaStep) > 0.00000001) {
		while (SSRdiff < 0.0 && (gammaTmp + gammaStep > 0)) {
			// update gamma
			gammaTmp += gammaStep;

			// fill x
			for (size_t p = 0; p < N; ++p) { X(p, 1) = exp(-gammaTmp * p); }
			betaHat = inv(X.t() * X) * X.t() * y;

			// calc sum of squares
			const arma::mat tmp = betaHat.t() * X.t() * y;
			const double SSRnew = sumYSquared - tmp(0, 0);
			SSRdiff             = SSRnew - SSRold;
			SSRold              = SSRnew;
		}
		// update alpha step
		gammaStep = -0.1 * gammaStep;
		SSRdiff   = -1.0;
	}
	// set parameters
	return {betaHat(0), betaHat(1), gammaTmp};
}

double calcLL(const std::vector<Probability> &Empiric, const std::array<double, 3> &Parameters) {
	double LL = 0.0;
	for (size_t p = 0; p < Empiric.size(); ++p) {
		const auto dExpMinusAlphaP = Parameters[1] * exp(-Parameters[2] * p);
		LL += Empiric[p].get() * log(Parameters[0] + dExpMinusAlphaP) +
			  (1. - Empiric[p].get()) * log(1.0 - Parameters[0] - dExpMinusAlphaP);
	}
	return LL;
}

auto fillFAndJacobian(const std::vector<Probability> &Empiric, const std::array<double, 3> &Parameters) {
	arma::vec::fixed<3> F;
	arma::mat::fixed<3,3> J;

	for (size_t p = 0; p < Empiric.size(); ++p) {
		// exp
		const auto expMinusAlphaP  = exp(-Parameters[2] * p);
		const auto dExpMinusAlphaP = Parameters[1] * expMinusAlphaP;

		// first term
		//----------
		const auto tmp1     = Parameters[0] + dExpMinusAlphaP;
		const auto weight1  = Empiric[p].get() / tmp1;
		const auto weightJ1 = weight1 / tmp1;

		// add to F
		F(0) += weight1;
		F(1) += weight1 * expMinusAlphaP;
		F(2) -= weight1 * p * dExpMinusAlphaP;

		// add to J -> only upper triangle, as it is symmetric
		J(0, 0) -= weightJ1;
		J(0, 1) -= weightJ1 * expMinusAlphaP;
		J(0, 2) += weightJ1 * p * dExpMinusAlphaP;
		J(1, 1) -= weightJ1 * expMinusAlphaP * expMinusAlphaP;
		J(1, 2) -= weightJ1 * p * Parameters[0] * expMinusAlphaP;
		J(2, 2) += weightJ1 * p * p * Parameters[0] * dExpMinusAlphaP;

		// second term
		//-----------
		const auto tmp2     = (1.0 - Parameters[0] - dExpMinusAlphaP);
		const auto weight2  = (1. - Empiric[p].get()) / tmp2;
		const auto weightJ2 = weight2 / tmp2;

		// add to F
		F(0) -= weight2;
		F(1) -= weight2 * expMinusAlphaP;
		F(2) += weight2 * p * dExpMinusAlphaP;

		// add to J -> only upper triangle, as it is symmetric
		J(0, 0) -= weightJ2;
		J(0, 1) -= weightJ2 * expMinusAlphaP;
		J(0, 2) += weightJ2 * p * dExpMinusAlphaP;
		J(1, 1) -= weightJ2 * expMinusAlphaP * expMinusAlphaP;
		J(1, 2) += weightJ2 * p * (1.0 - Parameters[0]) * expMinusAlphaP;
		J(2, 2) -= weightJ2 * p * p * (1.0 - Parameters[0]) * dExpMinusAlphaP;
	}

	// now fill in lower triangle of J
	J(1, 0) = J(0, 1);
	J(2, 0) = J(0, 2);
	J(2, 1) = J(1, 2);
	return std::make_tuple(F, J);
}

auto estimateWithNewtonRaphson(const std::vector<Probability> &Empiric, std::array<double, 3> Parameters) {
	const double epsilon = parameters().get<double>("epsilon", 0.001);
	logfile().list("Will consider the Newton-Raphson algorithm to have converged if the likelihood difference < " +
	               toString(epsilon) + ". (parameter 'epsilon')");
	const double numNRIterations = parameters().get<int>("numNR", 100);
	logfile().list("Will run up to " + toString(numNRIterations) + " Newton-Raphson iterations. (parameter 'numNR')");
	// Conduct Newton-Raphson to refine
	//----------------------------------
	double oldLL = impl::calcLL(Empiric, Parameters);
	for (size_t _ = 0; _ < numNRIterations; ++_) {
		auto [F, J] = impl::fillFAndJacobian(Empiric, Parameters);

		arma::vec::fixed<3> JxF;
		if (!solve(JxF, J, F)) {
			DEVERROR("Issue solving JxF in TPMDTable::fitExponentialModel!");
		}

		// estimate new params
		double lambda   = 1;
		double deltaLL  = 0;
		while (lambda > 1e-20) {
			std::array<double, 3> newParams;
			for (size_t i = 0; i < newParams.size(); ++i) newParams[i] = Parameters[i] - lambda * JxF(i);

			// calculate LL at new location
			const auto LL = impl::calcLL(Empiric, newParams);
			deltaLL = LL - oldLL;
			if (deltaLL > 0) {
				oldLL = LL;
				std::copy(newParams.begin(), newParams.end(), Parameters.begin());
				break;
			}

			// else backtrace
			lambda /= 2;
		}
		if (deltaLL < epsilon || lambda <= 1e20) break;
	}
	return std::make_tuple(Parameters[1], Parameters[2], Parameters[0]);
}

std::vector<Probability> makeEmpiric(const std::vector<double> &From_to, const std::vector<double> &To_from) {
	std::vector<Probability> values;
	values.reserve(From_to.size()); // include extra bin for sites beyond size (available in PMDTables)

	for (size_t p = 0; p < From_to.size(); ++p) {
		const double forward  = From_to[p]; // e.g. C -> T
		const double backward = To_from[p]; // e.g. T -> C

		if (forward == 0.) {
			values.push_back(0.);
		} else {
			values.push_back(std::max(0.0, (forward - backward) / (1.0 - backward)));
		}
	}
	return values;
}

} // namespace impl

//---------------------------------------------------------------
// TPMDFunctionNoPMD
//---------------------------------------------------------------
TNo::TNo(std::string_view string) {
	const auto params = impl::parseParameters<double>(string);
	if (params.size() != 0) {
		UERROR("Cannot initialize PMD function '", name, "': expected 0 but found ", params.size(), " parameters!");
	}
}

//---------------------------------------------------------------
// TPMDFunctionExponential
//--------------------------------------------------------------
TExponential::TExponential(std::string_view string) {
	constexpr size_t nParams = 4;
	const auto params = impl::parseParameters<double>(string);
	if (params.empty()) {
		// parameters missing: set to no PMD
		_a = _b = _c = 0.;
		_values.push_back(0.);
	} else {
		// parameters are provided
		if (params.size() != nParams) {
			UERROR("Cannot initialize PMD function '", name, "': expected", nParams, "(", example, ") but found ",
				   params.size(), " parameters!");
		}
		const auto N = std::lround(params[0]);
		_a           = params[1];
		_b           = params[2];
		_c           = params[3];

		if (N == 0)  UERROR("Cannot initialize PMD function '", name, "': last position must be > 0!"); 
		if (_a < 0.0) UERROR("Cannot initialize PMD function '", name, "': a must be > 0!");
		if (_b < 0.0) UERROR("Cannot initialize PMD function '", name, "': b must be > 0!");
		_fillPMDProbabilities(N + 1);
	}
}

void TExponential::_fillPMDProbabilities(size_t N) {
	_values.resize(N);
	for (size_t p = 0; p < N; ++p) { _values[p] = _a * exp(-_b * p) + _c; }
}


bool TExponential::learn(const std::vector<double> &From_to, const std::vector<double> &To_from) {
	logfile().list("Learning exponential pattern");

	_values.clear();

	if (From_to.size() < 3) {
		_a = 0.;
		_b = 0.;
		_c = 0.;
		_values.push_back(0.);
		return false;
	}

	try {
		// get initial estimates via OLS
		const auto empiric    = impl::makeEmpiric(From_to, To_from);
		const auto Parameters = impl::initialEstimatesOLS(empiric);
		std::tie(_a, _b, _c)  = impl::estimateWithNewtonRaphson(empiric, Parameters);
	} catch (...) {
		_a = 0.;
		_b = 0.;
		_c = 0.;
		_values.push_back(0.);
		return false;
	}

	if (_a < 0 || _b < 0) {
		_a = 0.;
		_b = 0.;
		_c = 0.;
		_values.push_back(0.);
		return false;
	}

	logfile().conclude(_a, "*exp(-", _b, "*p) + ", _c);

	_values.clear();
	_fillPMDProbabilities(From_to.size() + 1);

	// check if pattern is negativ
	if (_values.back() < 0) {
		_a = 0.;
		_b = 0.;
		_c = 0.;
		_values.push_back(0.);
		return false;
	}
	return true;
}

Probability TExponential::prob(uint16_t pos) const noexcept {
	// Note: distance is zero based!
	// model is fit up to _lastPosition. We assume constant PMD after that
	return pos < _values.size() ? _values[pos] : _values.back();
}

//---------------------------------------------------------------
// TPMDFunctionEmpiric
//---------------------------------------------------------------
TEmpiric::TEmpiric(std::string_view string) : _values(impl::parseParameters<Probability>(string)) {
	if (_values.empty()) {
		// parameters missing: set to no PMD
		_values.emplace_back(0.);
	} else {
		// parameters are provided
		for (auto &d : _values) {
			if (d < 0.0 || d > 1.0) {
				UERROR("Cannot initialize post mortem damage function '", name,
					   "': some probabilities are outside [0,1]!");
			}
		}
	}
}

bool TEmpiric::learn(const std::vector<double> &From_to, const std::vector<double> &To_from) {
	logfile().list("Learning empiric pattern");
	// resize parameters
	_values = impl::makeEmpiric(From_to, To_from);
	if (_values.empty()) {
		_values.emplace_back(0.);
		return false;
	}
	return true;
}

Probability TEmpiric::prob(uint16_t pos) const noexcept {
	return pos < _values.size() ? _values[pos] : _values.back();
}


TFunction* makeFunction(std::string_view pmdString) {
	// Options are
	//  none
	//  Empiric[0.5,0.3,...]
	//  Exponential[a,b,c]

	const auto name = coretools::str::readBefore(pmdString, '[');

	if (name == TNo::name) return new TNo(pmdString);
	if (name == TExponential::name) return new TExponential(pmdString);
	if (name == TEmpiric::name) return new TEmpiric(pmdString);

	UERROR("Cannot initialize PMD function: unknown function '", name, "'!. Use either ", TNo::name, ", ",
		   TExponential::name, " or ", TEmpiric::name + ".");
}

} // namespace GenotypeLikelihoods
