/*
 * TEpsilon.cpp
 *
 *  Created on: Jul 19, 2022
 *      Author: Andreas
 */

#include "TEpsilon.h"
#include "SequencingError/TCovariate.h"
#include "SequencingError/TFunction.h"
#include "stringFunctions.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

using namespace coretools::str;

namespace impl {

auto parseFunctions(const std::string &Defs) {
	struct Voldemort {
		std::string covariate;
		std::string function;
		Voldemort(std::string Covariate, std::string Function): covariate(std::move(Covariate)), function(std::move(Function)) {}
	};
	std::vector<Voldemort> functions;
	std::vector<std::string> defs;
	coretools::str::fillContainerFromString(Defs, defs, ';', true);

	for (const auto& def: defs) {
		const auto pos  = def.find(':');
		if (pos == std::string::npos) functions.emplace_back("", def);
		else functions.emplace_back(def.substr(0, pos), def.substr(pos + 1));
	}
	return functions;
}

auto parseFunction(const std::string &str) {
	constexpr auto format =
	    "Expected format is TYPE(ARGS)[BETAS], where (ARGS) is only required for some TYPE and [BETAS] is optional.";
	std::string type;
	std::vector<std::string> args;
	std::vector<std::string> betas;

	// split string into parameters and values
	if (const auto pos = str.find('['); pos != std::string::npos) {
		// extract type
		type = str.substr(0, pos);

		// extract parameters
		const auto pos2 = str.find(']');
		if (pos == std::string::npos) { throw "Wrong format for recal function '" + str + "': missing ']'! " + format; }
		coretools::str::fillContainerFromStringAny(str.substr(pos + 1, pos2 - pos - 1), betas, ",", true);
	} else {
		type = str;
	}

	// name might contain args
	if (const auto pos = type.find('('); pos != std::string::npos) {
		// extract parameters
		const auto pos2 = str.find(')');
		if (pos == std::string::npos) { throw "Wrong format for recal function '" + str + "': missing ')'! " + format; }
		coretools::str::fillContainerFromStringAny(type.substr(pos + 1, pos2 - pos - 1), args, ",", true);

		// extract type
		type = type.substr(0, pos);
	}
	return std::make_tuple(type, args, betas);
}

template<typename Covariate> TFunction *makeCovFunction(const std::string &Function, size_t FirstParameterIndex) {
	const auto [type, args, betas] = parseFunction(Function);
	if (type == TPolynomial<1, Covariate>::name) {
		size_t O = 0;
		if (betas.empty()) {
			if (args.size() != 1) throw "You must specify betas or order of Polynomial function!";
			fillFromString(args.front(), O);
		} else {
			if (((args.size() == 1) && (convertStringCheck<size_t>(args.front()) != betas.size()))
				|| !args.empty())
				throw "Number of betas does not correspond to order";
			O = betas.size();
		}
		if (O == 1) return new TPolynomial<1, Covariate>(FirstParameterIndex, betas);
		if (O == 2) return new TPolynomial<2, Covariate>(FirstParameterIndex, betas);
		if (O == 3) return new TPolynomial<3, Covariate>(FirstParameterIndex, betas);
		if (O == 4) return new TPolynomial<4, Covariate>(FirstParameterIndex, betas);
		if (O == 5) return new TPolynomial<5, Covariate>(FirstParameterIndex, betas);
		if (O == 6) return new TPolynomial<6, Covariate>(FirstParameterIndex, betas);
		throw "Only Polynomials up to order 6 can be used!";
	}
	if (!args.empty()) throw type + " cannot have arguments";
	if (type == TProbit<Covariate>::name) { return new TProbit<Covariate>(FirstParameterIndex, betas); }
	if (type == TEmpiric<Covariate>::name) {
		if (Covariate::isIndexed)
			return new TIndexedEmpiric<Covariate>(FirstParameterIndex, betas);
		else
			return new TEmpiric<Covariate>(FirstParameterIndex, betas);
	}
	if (type == TIndexedEmpiric<Covariate>::name) { return new TIndexedEmpiric<Covariate>(FirstParameterIndex, betas); }

	throw "Function '" + type + "' does not exist!";
}

TFunction *makeFunction(const std::string &Covariate, const std::string &Function, size_t FirstParameterIndex) {
	// No covariate
	if (Covariate.empty()) {
		const auto [type, args, betas] = parseFunction(Function);
		if (type != TIntercept::name) throw "You must specify a covariate";
		if (!args.empty()) throw "Intercept cannot have arguments";
		return new TIntercept(FirstParameterIndex, betas);
	}

	if (Covariate == TCovariate_quality::name)
		return makeCovFunction<TCovariate_quality>(Function, FirstParameterIndex);
	if (Covariate == TCovariate_position::name)
		return makeCovFunction<TCovariate_position>(Function, FirstParameterIndex);
	if (Covariate == TCovariate_context::name)
		return makeCovFunction<TCovariate_context>(Function, FirstParameterIndex);
	if (Covariate == TCovariate_fragmentLength::name)
		return makeCovFunction<TCovariate_fragmentLength>(Function, FirstParameterIndex);
	if (Covariate == TCovariate_mappingQuality::name)
		return makeCovFunction<TCovariate_mappingQuality>(Function, FirstParameterIndex);

	throw "Covariate '" + Covariate + "' does not exist!";
}

constexpr coretools::Probability calcEpsilon(double eta) noexcept {
	if (eta > 23.03) return coretools::Probability(0.9999999999);
	if (eta < -23.03) return coretools::Probability(0.0000000001);

	return coretools::logistic(eta);
}

} // namespace impl

TEpsilon::TEpsilon(const std::string &Def) {
	// create functions
	size_t _numParameters = 0;
	const auto modelDef = impl::parseFunctions(Def);
	for (const auto &cov : modelDef) {
		_functions.emplace_back(impl::makeFunction(cov.covariate, cov.function, _numParameters));

		// add new parameters
		_numParameters += _functions.back()->numParameters();
	}

	// prepare Newton-Raphson variables
	_Jacobian.resize(_numParameters, _numParameters);
	_F.resize(_numParameters);
	_JxF.resize(_numParameters, 1);
}

TEpsilon::~TEpsilon() = default;

void TEpsilon::checkOrInit(const RecalEstimatorTools::TRecalDataTable &DataTable) {
	// these may change, recalculate
	size_t _numParameters     = 0;

	for (auto &fn : _functions) {
		if (!fn->checkOrInitValueRange(DataTable)) {
			throw "Function " + fn->typeString() + " does not cover full range of data";
		}
		_numParameters += fn->numParameters();
	}

	// prepare Newton-Raphson variables
	_Jacobian.resize(_numParameters, _numParameters);
	_F.resize(_numParameters);
	_JxF.resize(_numParameters, 1);
}

coretools::Probability TEpsilon::calcErrorRate(const BAM::TSequencedBase &base) const noexcept {
	// eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function
	double eta = 0.;
	for (const auto& fn: _functions) eta += fn->getEta(base);
	return impl::calcEpsilon(eta);
}

coretools::Probability TEpsilon::_calcErrorRate(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
											   std::vector<T2ndDerivative> &der2) const noexcept {
	// eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function
	double eta = 0.;
	for (const auto &fn : _functions) eta += fn->getEta(base, der1, der2);
	return impl::calcEpsilon(eta);
}


void TEpsilon::solveJxF() {
	// Need to copy numbers to other triangle in Jacobian, as only upper triangle is filled when parsing sites
	_Jacobian = arma::symmatu(_Jacobian);

	// scale F and J by 1/#sites
	_Jacobian = _Jacobian / _numSitesAdded;
	_F        = _F / _numSitesAdded;
	_maxF     = std::max(_F.max(), -_F.min());
	if (!solve(_JxF, _Jacobian, _F))
		UERROR("Issue solving JxF! This may be due to a lack of data. Consider adding more sites. Jacobian: ",
		       _Jacobian);

	// automatically reset
	_Jacobian.zeros();
	_F.zeros();
	_numSitesAdded = 0;
	_oldQ          = _Q; // set Q to value before proposal
}
void TEpsilon::propose(double lambda) {
	size_t index = 0;
	for (const auto &fn : _functions) fn->proposeNewParameters(_JxF, index, lambda);

	_Q = 0; // not valid anymore
}

bool TEpsilon::acceptOrReject() {
	if (_Q > _oldQ) {
		_Q = 0; // reset for next iteration
		return true;
	}

	// else
	for (const auto &fn : _functions) fn->rejectProposedParameters();
	_Q = 0.; // not valid anymore

	return false;
}
void TEpsilon::adjust() {
	//for (const auto &fn : _functions) _intercept.intercept() += cov.function->adjustParametersPostEstimation();
}

std::string TEpsilon::getDefinition() const noexcept {
	std::string modelString = _functions.front()->modelString();
	for (size_t i = 1; i < _functions.size(); ++i) {
		const auto &fn = _functions[i];
		modelString += ";" + fn->modelString();
	}
	return modelString;
}

} // namespace SequencingError
} // namespace GenotypeLikelihoods
