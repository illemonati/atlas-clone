/*
 * TEpsilon.cpp
 *
 *  Created on: Jul 19, 2022
 *      Author: Andreas
 */

#include "TEpsilon.h"
#include "TSequencingErrorCovariate.h"
#include "TSequencingErrorCovariateFunction.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

namespace impl {


auto parseModuleString(const std::string &str) {
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

TCovariate *covariate(const std::string &type) {
	if (type == TCovariate::name) return nullptr;
	if (type == TCovariate_quality::name) return new TCovariate_quality;
	if (type == TCovariate_position::name) return new TCovariate_position;
	if (type == TCovariate_context::name) return new TCovariate_context;
	if (type == TCovariate_fragmentLength::name) return new TCovariate_fragmentLength;
	if (type == TCovariate_mappingQuality::name) return new TCovariate_mappingQuality;
	throw "Unknown recalibration covariate '" + type;
}

template<size_t O = 6>
TFunction *poly(size_t order, const size_t FirstParameterIndex, const std::vector<std::string> &betas) {
	if constexpr (O == 0)
		UERROR("Polynomial Order must be at least 1");
	else {
		if (order > O) UERROR("Polynomial Order cannot be higher than ", O);
		if (order == O)
			return new TPolynomial<O>(FirstParameterIndex, betas);
		else
			return poly<O - 1>(order, FirstParameterIndex, betas);
	}
}

TFunction *function(const std::string &functionString, const size_t FirstParameterIndex) {
		const auto [type, args, betas] = parseModuleString(functionString);
		// create function
		if (type == TPolynomial<1>::name) {
			if (betas.empty() && args.size() != 1) {
				throw "Wrong number of arguments for polynomial recal function '" + functionString +
					"': expect one argument (order).";
			}
			const size_t order = betas.empty() ? coretools::str::convertStringCheck<uint32_t>(args[0]) : betas.size();
			return poly(order, FirstParameterIndex, betas);
	    }
	    if (type == TSpecific::name) {
			return new TSpecific(FirstParameterIndex, betas);
		}
		if (type == TSpecificMap::name) {
		    return new TSpecificMap(FirstParameterIndex, betas);
	    }
	    throw "Recalibration function '" + type + "' not valid for covariate!";
    }


constexpr coretools::Probability calcEpsilon(double eta) noexcept {
	if (eta > 23.03) return coretools::Probability(0.9999999999);
	if (eta < -23.03) return coretools::Probability(0.0000000001);

	return coretools::logistic(eta);
}

std::string parseIntercept(const std::string &Def) {
	if (!(Def.rfind("intercept", 0) == 0)) throw "you must define an intercept!";

	const auto pos_1 = Def.find('[');
	const auto pos_2 = Def.find(']');
	if (pos_1 == std::string::npos) {
		if (pos_2 != std::string::npos) throw "Unable to understand model string '" + Def + "': missing '['";
		return "0.0";
	}
	if (pos_2 == std::string::npos) throw "Unable to understand model string '" + Def + "': missing ']'";
	return Def.substr(pos_1 + 1, pos_2 - pos_1 - 1);
}

auto parseCovariates(const std::string &Defs) {
	std::vector<std::pair<std::string, std::string>> covariates;
	std::vector<std::string> defs;
	coretools::str::fillContainerFromString(Defs, defs, ';', true);

	// skip first as intercept is already parsed
	for (size_t i = 1; i < defs.size(); ++i) {
		const auto &def = defs[i];
		const auto pos  = def.find(':');
		covariates.push_back({def.substr(0, pos), def.substr(pos + 1)});
	}
	return covariates;
}
} // namespace impl

TEpsilon::TEpsilon(const std::string &Def) : _intercept(0, impl::parseIntercept(Def)) {
	// create covariates
	_numParameters     = _intercept.numParameters();
	_num1stDerivatives = _intercept.numNonZeroFirstDerivatives();
	_num2ndDerivatives = _intercept.numNonZeroFirstDerivatives();

	const auto modelDef = impl::parseCovariates(Def);
	for (const auto &cov : modelDef) {
		// create function for each covariate
		if (cov.first == TCovariate::name) continue;

		_covariates.push_back(
			TCovariateModel{impl::covariate(cov.first), impl::function(cov.second, _numParameters)});

		// add new parameters
		_numParameters += _covariates.back().function->numParameters();
		_num1stDerivatives += _covariates.back().function->numNonZeroFirstDerivatives();
		_num2ndDerivatives += _covariates.back().function->numNonZeroSecondDerivatives();
	}

	// prepare Newton-Raphson variables
	_Jacobian.resize(_numParameters, _numParameters);
	_F.resize(_numParameters);
	_JxF.resize(_numParameters, 1);
}


void TEpsilon::checkOrInit(const RecalEstimatorTools::TRecalDataTable &DataTable) {
	// these may change, recalculate
	_numParameters     = _intercept.numParameters();
	_num1stDerivatives = _intercept.numNonZeroFirstDerivatives();
	_num2ndDerivatives = _intercept.numNonZeroFirstDerivatives();

	for (auto &cov : _covariates) {
		if (!cov.function->checkOrInitValueRange(cov.covariate->range(DataTable))) {
			throw "Function " + cov.function->typeString() + " does not cover full range of data of covariate " +
				   cov.covariate->typeString() + '\n';
		}

		// these may change, recalculate
		_numParameters     += _covariates.back().function->numParameters();
		_num1stDerivatives += _covariates.back().function->numNonZeroFirstDerivatives();
		_num2ndDerivatives += _covariates.back().function->numNonZeroSecondDerivatives();
	}
}

coretools::Probability TEpsilon::calcErrorRate(const BAM::TSequencedBase &base) const noexcept {
	// eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function
	auto eta = _intercept.getEtaTerm();
	for (const auto &cov : _covariates) eta += cov.function->getEtaTerm(cov.covariate->extract(base));

	return impl::calcEpsilon(eta);
}

void TEpsilon::addToEpsilon(const BAM::TSequencedBase &base, coretools::Probability P_g_I_d,
							coretools::Probability P_bbar_I_gd, bool update) {
	if (base == genometools::Base::N) return;

	// get error rate
	const double eps   = calcErrorRate(base);
	const double eps_c = 1. - eps;

	// add Q
	_Q += P_g_I_d.get() * (P_bbar_I_gd.get() * log(eps_c) + P_bbar_I_gd.complement().get() * log(eps));

	if (!update) return;

	// F and J
	const double w_ij = P_g_I_d * (eps_c - P_bbar_I_gd.get());

	static std::vector<T1stDerivative> der1st;
	der1st.clear();
	static std::vector<T2ndDerivative> der2nd;
	der2nd.clear();
	der1st.push_back(_intercept.get1stDerivatives());

	for (const auto &cov : _covariates) {
		for (size_t i = 0; i < cov.function->numNonZeroFirstDerivatives(); ++i) {
			der1st.push_back(cov.function->get1stDerivatives(cov.covariate->extract(base), i));
			for (size_t j = i; j < cov.function->numNonZeroSecondDerivatives(); ++j) {
				der2nd.push_back(cov.function->get2ndDerivatives(cov.covariate->extract(base), i, j));
			}
		}
	}

	// add first derivatives
	for (auto dm = der1st.begin(); dm != der1st.end(); ++dm) {
		// add to F
		_F(dm->index) += w_ij * dm->derivative;

		// add to J
		_Jacobian(dm->index, dm->index) -= eps_c * eps * dm->derivative * dm->derivative;
		for (auto dn = dm + 1; dn != der1st.end(); ++dn) {
			_Jacobian(dm->index, dn->index) -= eps_c * eps * dm->derivative * dn->derivative;
		}
	}

	// add second derivatives to Jacobian (happens to have the same weigth as F!)
	for (auto &dmn : der2nd) _Jacobian(dmn.index1, dmn.index2) += w_ij * dmn.derivative;

	++_numSitesAdded;
}

void TEpsilon::solveJxF() {
	// Need to copy numbers to other triangle in Jacobian, as only upper triangle is filled when parsing sites
	for (size_t i = 0; i < _numParameters - 1; ++i) {
		for (size_t j = i + 1; j < _numParameters; ++j) {
			// copy from upper triangle to lower triangle
			_Jacobian(j, i) = _Jacobian(i, j);
		}
	}

	// scale F and J by 1/#sites
	_Jacobian = _Jacobian / (double)_numSitesAdded;
	_F        = _F / (double)_numSitesAdded;
	_maxF     = std::max(_F.max(), -_F.min());
	if (!solve(_JxF, _Jacobian, _F))
		UERROR("Issue solving JxF! This may be due to a lack of data. Consider adding more sites. Jacobian: ",
			   _Jacobian);

	// automatically reset
	_Jacobian.zeros();
	_F.zeros();
	_numSitesAdded = 0;
}
void TEpsilon::proposeNewParameters(double lambda) {
	uint16_t index = 0;
	_intercept.proposeNewParameters(_JxF, index, lambda);
	for (const auto &cov : _covariates) cov.function->proposeNewParameters(_JxF, index, lambda);
}
bool TEpsilon::acceptProposedParametersBasedOnQ() {
	if (_Q > _oldQ) { return true; }

	// else
	_intercept.rejectProposedParameters();
	for (const auto &cov : _covariates) cov.function->rejectProposedParameters();
	_Q = _oldQ;

	return false;
}
void TEpsilon::adjustParametersPostEstimation() {
	for (const auto &cov : _covariates) _intercept.intercept() += cov.function->adjustParametersPostEstimation();
}

std::string TEpsilon::getDefinition() const noexcept {
	std::string modelString = _intercept.modelString();
	for (const auto &cov : _covariates) {
		modelString += ";" + cov.covariate->typeString() + ":" + cov.function->modelString();
	}
	return modelString;
}

} // namespace SequencingError
} // namespace GenotypeLikelihoods
