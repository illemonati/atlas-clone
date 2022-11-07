/*
 * TEpsilon.cpp
 *
 *  Created on: Jul 19, 2022
 *      Author: Andreas
 */

#include "TEpsilon.h"
#include "SequencingError/TCovariate.h"
#include "SequencingError/TFunction.h"
#include "coretools/Main/TError.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/splitters.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Strings/toString.h"
#include <algorithm>
#include <iterator>

namespace GenotypeLikelihoods {
namespace SequencingError {

using namespace coretools::str;

namespace impl {

auto parseFunctions(std::string_view Defs) {
	struct Voldemort {
		std::string_view covariate;
		std::string_view function;
		Voldemort(std::string_view Covariate, std::string_view Function): covariate(Covariate), function(Function) {}
	};
	std::vector<Voldemort> functions;

	for (auto def: TSplitter(Defs, ';')) {
		const auto pos  = def.find(':');
		if (pos == std::string_view::npos) functions.emplace_back("", def);
		else functions.emplace_back(def.substr(0, pos), def.substr(pos + 1));
	}
	return functions;
}

auto parseFunction(std::string_view str) {
	const auto beg = str.find('[');
	if (beg == std::string_view::npos) { // default arguments
		return std::make_pair(str, TSplitter<>("", ','));
	}
	const auto end = str.find(']', beg);
	if (end == std::string_view::npos) {
		UERROR("Wrong format for recal function '", str, "': missing ']'! ",
			   "Expected format is TYPE[BETAS], where [BETAS] is optional.");
	}
	return std::make_pair(str.substr(0, beg), TSplitter<>(str.substr(beg + 1, end - beg - 1), ','));
}

template<typename Covariate> TFunction *makeCovFunction(std::string_view Function, size_t FirstParameterIndex) {
	auto [type, Spl] = parseFunction(Function);

	if (stringStartsWith(type, TPolynomial<1, Covariate>::name)) {
		if (type.size() > TPolynomial<1, Covariate>::name.size() + 1) UERROR("Unknow function name: ", type, "!");
		std::array<double, 9> bs{};
		size_t o = 0;
		while(!Spl.empty()) {
			fromString<true>(Spl.front(), bs[o]);
			Spl.popFront();
			++o;
		}
		if (type.size() == TPolynomial<1, Covariate>::name.size() + 1) {
			size_t po = type.back() - '0';
			if (o == 0) o = po; // no values given
			else if (o != po) UERROR("Defined polynomial function of order ", po, " but gave ", o, " parameters!");
		}
		TFunction* fn;
		switch (o) {
		case 1: fn = new TPolynomial<1, Covariate>(FirstParameterIndex); break;
		case 2: fn = new TPolynomial<2, Covariate>(FirstParameterIndex); break;
		case 3: fn = new TPolynomial<3, Covariate>(FirstParameterIndex); break;
		case 4: fn = new TPolynomial<4, Covariate>(FirstParameterIndex); break;
		case 5: fn = new TPolynomial<5, Covariate>(FirstParameterIndex); break;
		case 6: fn = new TPolynomial<6, Covariate>(FirstParameterIndex); break;
		case 7: fn = new TPolynomial<7, Covariate>(FirstParameterIndex); break;
		case 8: fn = new TPolynomial<8, Covariate>(FirstParameterIndex); break;
		case 9: fn = new TPolynomial<9, Covariate>(FirstParameterIndex); break;
		default: UERROR("Only Polynomials from order 1 to 9 can be used!");
		}
		size_t i = 0;
		for (auto &beta : *fn) {
			beta = bs[i];
			++i;
		}
		return fn;
	}
	/*
	if (type == TProbit<Covariate>::name) {
		return new TProbit<Covariate>(FirstParameterIndex, Spl);
	}
	if (type == TEmpiric<Covariate>::name) {
		if (Covariate::isIndexed)
			return new TIndexedEmpiric<Covariate>(FirstParameterIndex, Spl);
		else
			return new TEmpiric<Covariate>(FirstParameterIndex, Spl);
	}
	*/

	UERROR("Function '", type, "' does not exist!");
}

TFunction *makeFunction(std::string_view Covariate, std::string_view Function, size_t FirstParameterIndex) {
	// No covariate
	if (Covariate.empty()) {
		const auto [type, betas] = parseFunction(Function);
		if (type != TIntercept::name) UERROR("You must specify a covariate");
		TFunction *in = new TIntercept(FirstParameterIndex);
		if (betas.empty()) return in;
		fromString<true>(betas.front(), *in->begin());
		return in;
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

	UERROR("Covariate '", Covariate, "' does not exist!");
}

constexpr coretools::Probability calcEpsilon(double eta) noexcept {
	if (eta > 23.03) return coretools::Probability(0.9999999999);
	if (eta < -23.03) return coretools::Probability(0.0000000001);

	return coretools::logistic(eta);
}

} // namespace impl

TEpsilon::TEpsilon(std::string_view Def) {
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
	_oldBetas.clear();
	for (const auto &fn : _functions) {
		for (auto& beta: *fn) {
			_oldBetas.push_back(beta);
			beta -= lambda * _JxF(index++);
		}
	}

	_Q = 0; // not valid anymore
}

bool TEpsilon::acceptOrReject() {
	if (_Q > _oldQ) {
		_Q = 0; // reset for next iteration
		return true;
	}

	// else
	auto old = _oldBetas.begin();
	for (const auto &fn : _functions) {
		for (auto& beta: *fn) {
			beta = *old;
			++old;
		}
	}
	_Q = 0.; // not valid anymore

	return false;
}
void TEpsilon::adjust() {
	//for (const auto &fn : _functions) _intercept.intercept() += cov.function->adjustParametersPostEstimation();
}

std::string TEpsilon::getDefinition() const noexcept {
	/*std::string modelString = _functions.front()->modelString();
	for (size_t i = 1; i < _functions.size(); ++i) {
		const auto &fn = _functions[i];
		modelString.push_back(';');
		modelString.append(fn->modelString());
	}
	return modelString;*/
	return std::accumulate(_functions.begin() + 1, _functions.end(), _functions.front()->modelString(),
						   [](auto tot, auto& f) { return tot + ";" + f->modelString(); });
}

} // namespace SequencingError
} // namespace GenotypeLikelihoods
