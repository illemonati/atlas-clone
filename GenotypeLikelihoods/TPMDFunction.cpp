#include "TPMDFunction.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"

namespace GenotypeLikelihoods {

using namespace coretools::str;
using coretools::Probability;
using coretools::instances::logfile;
using coretools::instances::parameters;
using genometools::Base;

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
}

//---------------------------------------------------------------
// TPMDFunctionNoPMD
//---------------------------------------------------------------
TPMDFunctionNoPMD::TPMDFunctionNoPMD(std::string_view string) {
	const auto params = impl::parseParameters<double>(string);
	if (params.size() != 0) {
		UERROR("Cannot initialize PMD function '", name, "': expected 0 but found ", params.size(), " parameters!");
	}
}

//---------------------------------------------------------------
// TPMDFunctionExponential
//--------------------------------------------------------------
TPMDFunctionExponential::TPMDFunctionExponential(std::string_view string) {
	constexpr size_t nParams = 4;
	const auto params = impl::parseParameters<double>(string);
	if (params.empty()) {
		// parameters missing: set to no PMD
		_lastPosition = 0;
		_a = _b = _c = 0.;
	} else {
		// parameters are provided
		if (params.size() != nParams) {
			UERROR("Cannot initialize PMD function '", name, "': expected", nParams, "(", example, ") but found ",
				   params.size(), " parameters!");
		}
		_lastPosition = std::lround(params[0]);
		_a = params[1];
		_b = params[2];
		_c = params[3];

		if (_lastPosition == 0)  UERROR("Cannot initialize PMD function '", name, "': last position must be > 0!"); 
		if (_b < 0.0) UERROR("Cannot initialize PMD function '", name, "': b must be > 0!");
		if (_c < 0.0) UERROR("Cannot initialize PMD function '", name, "': c must be > 0!");
	}
	_fillPMDProbabilities();
}

void TPMDFunctionExponential::_fillPMDProbabilities() {
	_values.resize(_lastPosition + 1);
	for (size_t p = 0; p < _values.size(); ++p) {
		_values[p] = _a * exp(-_b * p) + _c;
	}
}

	void TPMDFunctionExponential::parseEstimationParameters(TPMDEstimationParameters &EstimationParameters) {
	if (EstimationParameters.find(epsilon) == EstimationParameters.end()) {
		double eps = parameters().getParameterWithDefault<double>(epsilon, 0.001);
		EstimationParameters.emplace(epsilon, eps);
		logfile().list("Will consider the Newton-Raphson algorithm to have converged if the likelihood difference < " +
			      toString(eps) + ". (parameter '" + epsilon + "')");
	}

	if (EstimationParameters.find(numNR) == EstimationParameters.end()) {
		double numNRIterations = parameters().getParameterWithDefault<int>(numNR, 100);
		EstimationParameters.emplace(numNR, numNRIterations);
		logfile().list("Will run up to " + toString(numNRIterations) + " Newton-Raphson iterations. (parameter '" +
			      numNR + "')");
	}
}

void TPMDFunctionExponential::_initialEstimatesOLS(const countVec &pmdCounts, const countVec &pmdSums,
						   std::array<double, 3> &Parameters) {
	// fill vector y to fit using OLS
	const size_t N = _lastPosition + 1;
	arma::vec y(N);
	double sumYSquared = 0.0;
	for (int p = 0; p <= _lastPosition; ++p) {
		y(p) = (double)pmdCounts[p]/pmdSums[p];
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
		while (SSRdiff < 0.0) {
			// update gamma
			gammaTmp += gammaStep;

			// fill x
			for (size_t p = 0; p < N; ++p) {
				X(p, 1) = exp(-gammaTmp*p);
			}
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
	Parameters = {betaHat(0), betaHat(1), gammaTmp};
}

void TPMDFunctionExponential::_fillFAndJacobian(arma::vec &F, arma::mat &J, const countVec &pmdCounts,
						const countVec &pmdSums, const std::array<double, 3> &Parameters) {
	F.zeros();
	J.zeros();

	for (int p = 0; p <= _lastPosition; ++p) {
		// exp
		const auto expMinusAlphaP  = exp(-Parameters[2] * p);
		const auto dExpMinusAlphaP = Parameters[1] * expMinusAlphaP;

		// first term
		//----------
		const auto tmp1     = Parameters[0] + dExpMinusAlphaP;
		const auto weight1  = pmdCounts[p] / tmp1;
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
		const auto weight2  = (pmdSums[p] - pmdCounts[p]) / tmp2;
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
}

double TPMDFunctionExponential::_calcLL(const countVec &pmdCounts, const countVec &pmdSums,
					const std::array<double, 3> &Parameters) {
	double LL = 0.0;
	for (int p = 0; p <= _lastPosition; ++p) {
		const auto dExpMinusAlphaP = Parameters[1] * exp(-Parameters[2] * p);
		LL += pmdCounts[p] * log(Parameters[0] + dExpMinusAlphaP) +
		      (pmdSums[p] - pmdCounts[p]) * log(1.0 - Parameters[0] - dExpMinusAlphaP);
	}
	return LL;
}

void TPMDFunctionExponential::_estimateWithNewtonRaphson(const countVec &pmdCounts, const countVec &pmdSums,
							 std::array<double, 3> &Parameters, uint32_t numNRIterations,
							 double epsilon) {
	// Conduct Newton-Raphson to refine
	//----------------------------------
	double oldLL = _calcLL(pmdCounts, pmdSums, Parameters);

	for (size_t _ = 0; _ < numNRIterations; ++_) {
		arma::mat J(3, 3);
		arma::vec F(3);
		arma::mat JxF;
		_fillFAndJacobian(F, J, pmdCounts, pmdSums, Parameters);

		if (!solve(JxF, J, F)) {
			DEVERROR("Issue solving JxF in TPMDTable::fitExponentialModel!");
		}

		// estimate new params
		std::array<double, 3> newParams;
		for (size_t i = 0; i < newParams.size(); ++i) newParams[i] = Parameters[i] - JxF(i);

		// calculate LL at new location
		const auto LL = _calcLL(pmdCounts, pmdSums, newParams);

		// check if we accept or backtrack
		if (LL > oldLL) {
			oldLL = LL;
			// store new params
			for (size_t i = 0; i < Parameters.size(); ++i) Parameters[i] = newParams[i];

			// check if we stop NR
			if (LL - oldLL < epsilon) break;
		}
	}
}

void TPMDFunctionExponential::learn(const TPMDTable &Table, const Base &from, const Base &to,
				    const TPMDEstimationParameters &EstimationParameters) {
	logfile().list("Learning exponential pattern");
	// extract counts in PMD direction and the inverse direction
	const countVec &pmdCounts = Table[from][to];
	const countVec &pmdSums   = Table.sums(from);
	const countVec &invCounts = Table[to][from];
	const countVec &invSums   = Table.sums(from);

	// find last entry with counts
	_lastPosition = pmdCounts.size() - 1;
	while (_lastPosition >= 9 && pmdCounts[_lastPosition] < 100) _lastPosition --;

	// Check if we have sufficient data
	if (_lastPosition <= 9)
		UERROR("Not sufficient data to fit exponential PMD model: less than the ten first positions have > 100 data "
			   "points!\nConsider pooling read groups (parameter poolReadGroups).");

	const auto last = pmdSums.cbegin() + _lastPosition;
	if (std::find(pmdSums.cbegin(), last, 0) != last)
		UERROR("Not sufficient data to fit exponential PMD model: no observations for some reference "
			   "alleles!<nConsider reducing the relevant length (parameter length).");

	// get initial estimates via OLS
	std::array<double, 3> Parameters;
	_initialEstimatesOLS(pmdCounts, pmdSums, Parameters);


	// run Newton-Raphson
	_estimateWithNewtonRaphson(pmdCounts, pmdSums, Parameters,
				   EstimationParameters.at(numNR),
				   EstimationParameters.at(epsilon));


	// transform parameters
	// the exponential PMD model is f(C->T) = mu + (1-mu) *[ a*exp(-b * position) + c ]
	// but we fitted f(C->T) = alpha + beta * exp(-gamma * position).
	// Hence we have:
	//  mu is estimated from T->C transitions
	//  a =  beta / (1 - mu)
	//  b = gamma
	//  c = (alpha - mu) / (1 - mu)
	const double mu = std::accumulate(invCounts.cbegin(), invCounts.cbegin() + _lastPosition + 1, 0.) /
			  std::accumulate(invSums.cbegin(), invSums.cbegin() + _lastPosition + 1, 0);

	// store parameters, including lastPosition
	_a = Parameters[1] / (1.0 - mu);
	_b = Parameters[2];
	_c = (Parameters[0] - mu) / (1.0 - mu);

	logfile().conclude(_a, "*exp(-", _b, "*p) + ", _c);

	if (Parameters[1] < 0) {
		UERROR("Estimation resulted in a = ", _a,
			   " < 0! This is likely due to limited data. Consider pooling read groups (parameter poolReadGroups).");
	}

	_fillPMDProbabilities();


	// check if pattern is negativ
	if (_values[_lastPosition] < 0) {
		UERROR("Estimation resulted in negative PMD = ", _values[_lastPosition],
			   " at high positions!\nThis is likely be due to limited data. Consider pooling read groups (parameter "
			   "poolReadGroups).");
	}
}

Probability TPMDFunctionExponential::prob(uint16_t pos) const noexcept {
	// Note: distance is zero based!
	// model is fit up to _lastPosition. We assume constant PMD after that
	return pos < _lastPosition ? _values[pos] : _values[_lastPosition];
}

//---------------------------------------------------------------
// TPMDFunctionEmpiric
//---------------------------------------------------------------
	TPMDFunctionEmpiric::TPMDFunctionEmpiric(std::string_view string) : _values(impl::parseParameters<Probability>(string)) {
	if (_values.empty()) {
		// parameters missing: set to no PMD
		_values = {0.0};
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

void TPMDFunctionEmpiric::learn(const TPMDTable &Table, const Base &from, const Base &to,
				const TPMDEstimationParameters &) {
	logfile().list("Learning empiric pattern");
	// resize parameters
	_values.resize(Table.size()); // include extra bin for sites beyond size (available in PMDTables)

	// extract counts in PMD direction and the inverse direction
	const countVec &forwardCounts  = Table[from][to]; //e.g. C -> T
	const countVec &forwardSums    = Table.sums(from);
	const countVec &backwardCounts = Table[to][from]; //e.g. T -> C
	const countVec &backwardSums   = Table.sums(to);

	for (size_t p = 0; p < _values.size(); ++p) {
		if (forwardSums[p] == 0 || backwardSums[p] == 0) {
			_values[p] = 0.0;
		} else {
			double forward  = (double)forwardCounts[p] / forwardSums[p]; //e.g. C -> T
			double backward = (double)backwardCounts[p] / backwardSums[p]; //e.g. T -> C

			//forward = mu_CT + (1 - mu_CT) * PMD; mu_CT = mu_TC = backward
			_values[p] = std::max(0.0, (forward - backward) / (1.0 - backward));
		}
	}
}

Probability TPMDFunctionEmpiric::prob(uint16_t pos) const noexcept {
	return pos < _values.size() ? _values[pos] : _values.back();
}


TPMDFunction* makeFunction(std::string_view pmdString) {
	// Options are
	//  none
	//  Empiric[0.5,0.3,...]
	//  Exponential[a,b,c]

	const auto name = coretools::str::readBefore(pmdString, '[');

	if (name == TPMDFunctionNoPMD::name) return new TPMDFunctionNoPMD(pmdString);
	if (name == TPMDFunctionExponential::name) return new TPMDFunctionExponential(pmdString);
	if (name == TPMDFunctionEmpiric::name) return new TPMDFunctionEmpiric(pmdString);

	UERROR("Cannot initialize PMD function: unknown function '", name, "'!. Use either ", TPMDFunctionNoPMD::name, ", ",
		   TPMDFunctionExponential::name, " or ", TPMDFunctionEmpiric::name + ".");
}

} // namespace GenotypeLikelihoods
