/*
 * TRecalibrationEMModel.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TSequencingErrorModel.h"

#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <memory>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <utility>

#include "RecalEstimatorTools.h"
#include "TError.h"
#include "TGenotypeData.h"
#include "TRandomGenerator.h"
#include "TSequencedBase.h"
#include "TSequencingErrorCovariate.h"
#include "devtools.h"
#include "mathFunctions.h"
#include "probability.h"
#include "stringFunctions.h"
#include "toString.h"
#include "weakTypes.h"

namespace GenotypeLikelihoods {
namespace SequencingError {
using coretools::Probability;
using coretools::instances::randomGenerator;
using genometools::Base;
using namespace std::literals;

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

} // namespace impl

//*********************************************************
// TRecalibrationEMModelCovariateDefinition
// class to store model definition. Used when parsing files
//*********************************************************

TModelDefinition::TModelDefinition(const std::string &modelString, const std::string &rhoString) : rho(rhoString) {
	// split string
	std::vector<std::string> tmp;
	coretools::str::fillContainerFromString(modelString, tmp, ';', true);
	// loop over entries
	for (std::string s : tmp) {
		const auto pos = s.find(':');
		if (pos == std::string::npos) {
			// intercept?
			const auto pos_1 = s.find('[');
			const auto pos_2 = s.find(']');
			if (s.substr(0, pos_1) != "intercept")
				throw "Unable to understand model string '" + modelString +
					"': expecting an ':' or the function 'intercept'!";
			if (pos_1 == std::string::npos) {
				if (pos_2 != std::string::npos)
					throw "Unable to understand model string '" + modelString + "': missing '['";
				intercept = "0.0";
			} else {
				if (pos_2 == std::string::npos)
					throw "Unable to understand model string '" + modelString + "': missing ']'";
				intercept = s.substr(pos_1 + 1, pos_2 - pos_1 - 1);
			}
		} else {
			covariates.push_back({s.substr(0, pos), s.substr(pos + 1)});
		}
	}
}

//*********************************************************
// TRho
//*********************************************************

TRho::TRho(const std::string &def) {
	using coretools::str::toString;
	//"default" implies default rho
	if (def == "default") {
		return;
	}

	// otherwise: full matrix is provided
	std::vector<std::string> vec;
	std::string s = def;
	coretools::str::fillContainerFromString(def, vec, ';');
	if (vec.size() != 4) throw "Rho matrix has " + toString(vec.size()) + " instead of 4 rows!";

	// parse rows
	std::vector<double> r;
	for (Base a = Base::min; a < Base::max; ++a) {
		std::string &row = vec[index(a)];
		coretools::str::trimString(row, "()");
		coretools::str::fillContainerFromString(row, r, ',');
		if (r.size() != 4)
			throw "Rho matrix has " + toString(r.size()) + " instead of 4 columns for row " + toString(index(a) + 1) + "!";

		// fill
		for (Base b = Base::min; b < Base::max; ++b) {
			if (a == b) {
				rho[a][b] = 0.0;
			} else {
				rho[a][b] = r[index(b)];
			}
		}
	}
}

std::string TRho::getDefinition() const noexcept {
	using coretools::str::toString;
	return "-,"s + toString(rho[Base::A][Base::C]) + ',' + toString(rho[Base::A][Base::G]) + ',' + toString(rho[Base::A][Base::T]) + ';'
		+ toString(rho[Base::C][Base::A]) + ",-," + toString(rho[Base::C][Base::G]) + toString(rho[Base::C][Base::T]) + ';'
		+ toString(rho[Base::G][Base::A]) + ',' + toString(rho[Base::G][Base::C]) + ",-," + toString(rho[Base::G][Base::T]) + ';'
		+ toString(rho[Base::T][Base::A]) + ',' + toString(rho[Base::T][Base::C]) + ',' + toString(rho[Base::T][Base::G]) + ",-";
}

void TRho::add(genometools::Base base, coretools::Probability P_g_I_d, const TBaseLikelihoods &P_bbar_I_d) noexcept {
	for (auto b = Base::min; b < Base::max; ++b) {
		rho[b][base] += P_g_I_d*P_bbar_I_d[b];
	}
}

void TRho::estimate() noexcept {
	OUT(rho);
	for (Base a = Base::min; a < Base::max; ++a) {
		rho[a][a] = 0.0;
		double d  = 0.;
		for (const auto r : rho[a]) d += r;
		for (auto &r : rho[a]) r /= d;
	}
	OUT(rho);
}

//*********************************************************
// TModelNoRecal
//*********************************************************
Probability TModelNoRecal::getErrorRate(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) return Probability(1.0);
	return (Probability)base.originalQuality_phredInt;
}

genometools::PhredIntProbability TModelNoRecal::getPhredInt(const BAM::TSequencedBase &base) const noexcept {
	// Todo: change to maxProb() one available.
	if (base == Base::N) return genometools::PhredIntProbability(0);
	return base.originalQuality_phredInt;
}

TBaseLikelihoods TModelNoRecal::getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) {
		return TBaseLikelihoods{1.};
	}
	const auto eps = static_cast<Probability>(base.originalQuality_phredInt);
	TBaseLikelihoods baseLikelihoods{(1./3)*eps};
	baseLikelihoods[base.base] = eps.complement();
	return baseLikelihoods;
}

void TModelNoRecal::simulate(BAM::TSequencedBase &base) const noexcept {
	if (base.base == Base::N) return;

	const auto e = static_cast<Probability>(base.originalQuality_phredInt);
	if (randomGenerator().getRand() < e) {
		const int i = randomGenerator().getRand(1, 4); // 3 bases to choose from
		base.base   = Base((index(base.base) + i) % 4);
	}
}

//*********************************************************
// TModelRecal
//*********************************************************


TModelRecal::TModelRecal(const TModelDefinition &modelDef) : _rho(modelDef.rho), _intercept(0, modelDef.intercept) {
	// create covariates
	_functions.push_back(&_intercept);

	_numParameters = _intercept.numParameters();
	for (const auto & cov : modelDef.covariates) {
		// create function for each covariate
		if (cov.covariate == TCovariate::name) continue;

		// Andreas
		// auto cov = covriate(cov.covariate);
		// auto fn = cov.getFunction(cov.function, _numParameters);
		// Andreas
		_covariates.push_back(TCovariateModel{impl::covariate(cov.covariate), impl::function(cov.function, _numParameters)});
		_functions.push_back(_covariates.back().function.get());

		// add new parameters
		_numParameters += _covariates.back().function->numParameters();
	}

	// prepare Newton-Raphson variables
	_Jacobian.resize(numParameters(), numParameters());
	_F.resize(numParameters());
	_JxF.resize(numParameters(), 1);
}

TModelDefinition TModelRecal::getModelDefinition() const {
	return TModelDefinition(getCovariateDefinition(), getRhoDefinition());
}

std::string TModelRecal::getCovariateDefinition() const noexcept {
	std::string modelString = _intercept.modelString();
	for (const auto &cov : _covariates) {
		modelString += ";" + cov.covariate->typeString() + ":" +  cov.function->modelString();
	}
	return modelString;
}

//-------------------------------------------------
// functions to calculate error rates
//-------------------------------------------------

constexpr Probability calcEpsilon(double eta) noexcept {
	if (eta > 23.03) return Probability(0.9999999999);
	if (eta < -23.03) return Probability(0.0000000001);

	return coretools::logistic(eta);
}

Probability TModelRecal::_calcErrorRate(const BAM::TSequencedBase &base) const noexcept {
	// eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function
	auto eta = _intercept.getEtaTerm();
	for (const auto &cov : _covariates) eta += cov.function->getEtaTerm(cov.covariate->extract(base));

	return calcEpsilon(eta);
}

Probability TModelRecal::getErrorRate(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) return Probability::highest();
	return _calcErrorRate(base);
}

genometools::PhredIntProbability TModelRecal::getPhredInt(const BAM::TSequencedBase &base) const noexcept {
	// Todo: change to maxProb() one available.
	if (base == Base::N) return genometools::PhredIntProbability(0);
	return genometools::PhredIntProbability(_calcErrorRate(base));
}

TBaseLikelihoods TModelRecal::getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept {
	if (base == Base::N) {
		return TBaseLikelihoods{1.};
	}
	const auto e = _calcErrorRate(base);
	TBaseLikelihoods baseLikelihoods;
	for (auto b = Base::min; b < Base::max; ++b) baseLikelihoods[b] = e * _rho(b, base.base);
	baseLikelihoods[base.base] = e.complement();
	return baseLikelihoods;
}

void TModelRecal::simulate(BAM::TSequencedBase &base) const noexcept {
	if (base.base == Base::N) return;

	const auto e = _calcErrorRate(base);
	if (randomGenerator().getRand() < e) {
		const double r = randomGenerator().getRand();
		double cumul   = 0.;
		for (auto b = Base::min; b < Base::max; ++b) {
			cumul += _rho(b, base.base); //_rho(base.base, base.base) = 0
			if (r < cumul) {
				base.base = b;
				return;
			}
		}
	}
}

//-------------------------------------------------
// functions to estimate rho
//-------------------------------------------------
void TModelRecal::resetRho() noexcept { _rho.reset(); }

void TModelRecal::addToRho(const BAM::TSequencedBase &data, coretools::Probability P_g_I_d, const TBaseLikelihoods &P_bbar_I_d) noexcept {
	_rho.add(data.base, P_g_I_d, P_bbar_I_d);
}

void TModelRecal::estimateRho() noexcept { _rho.estimate(); }

//-------------------------------------------------
// functions for estimation
//-------------------------------------------------
void TModelRecal::checkOrInit(const RecalEstimatorTools::TRecalDataTable &DataTable) const {
	for (auto &cov : _covariates) {
		if (!cov.function->checkOrInitValueRange(cov.covariate->range(DataTable))) {
			throw "Function " + cov.function->typeString() + " does not cover full range of data of covariate " +
				   cov.covariate->typeString() + '\n';
		}
	}
}

void TModelRecal::resetQFJ() noexcept {
	_Jacobian.zeros();
	_F.zeros();
	_numSitesAdded  = 0;
	_NRStepAccepted = false;
	_oldQ = _Q;
	_Q    = 0.0;
}

void TModelRecal::addToQFJ(const BAM::TSequencedBase &base, Probability P_g_I_d, Probability P_bbar_I_gd) {
	// get error rate
	const double eps   = getErrorRate(base);
	const double eps_c = 1. - eps;

	// add Q
	_Q += P_g_I_d * (P_bbar_I_gd * eps_c + P_bbar_I_gd.complement() * eps);

	// F and J
	const double w_ij = P_g_I_d * (eps_c - P_bbar_I_gd);

	std::vector<T1stDerivative> der1st;
	std::vector<T2ndDerivative> der2nd;
	der1st.push_back(_intercept.get1stDerivatives());

	for (const auto & cov: _covariates) {
		for (size_t i = 0; i < cov.function->numNonZeroFirstDerivatives(); ++i) {
			der1st.push_back(cov.function->get1stDerivatives(cov.covariate->extract(base), i));
			for (size_t j = i+1; j < cov.function->numNonZeroSecondDerivatives(); ++j) {
				der2nd.push_back(cov.function->get2ndDerivatives(cov.covariate->extract(base), i, j));
			}
		}
	}

	// add first derivatives
	for (auto dm = der1st.begin(); dm != der1st.end(); ++dm) {
		// add to F
		_F(dm->index) += w_ij * dm->derivative;

		// add to J
		for (auto dn = dm; dn != der1st.end(); ++dn) {
			_Jacobian(dm->index, dn->index) -= eps_c * eps * dm->derivative * dn->derivative;
		}
	}

	// add second derivatives to Jacobian (happens to have the same weigth as F!)
	for (auto &dmn : der2nd) _Jacobian(dmn.index1, dmn.index2) += w_ij * dmn.derivative;

	++_numSitesAdded;
}


bool TModelRecal::solveJxF() {
	// Need to copy numbers to other triangle in Jacobian, as only upper triangle is filled when parsing sites
	for (int i = 0; i < (numParameters() - 1); ++i) {
		for (unsigned int j = i + 1; j < numParameters(); ++j) {
			// copy from upper triangle to lower triangle
			_Jacobian(j, i) = _Jacobian(i, j);
		}
	}

	// scale F and J by 1/#sites
	_Jacobian = _Jacobian / (double)_numSitesAdded;
	_F        = _F / (double)_numSitesAdded;

	OUT(_Jacobian);
	OUT(_F);

	// now solve J^-1 x F
	return solve(_JxF, _Jacobian, _F);
}

void TModelRecal::proposeNewParameters(double lambda) {
	if (!_NRStepAccepted) {
		uint16_t index = 0;
		for (const auto f : _functions) { f->proposeNewParameters(_JxF, index, lambda); }
	}
}

bool TModelRecal::acceptProposedParametersBasedOnQ() {
	OUT(_Q);
	OUT(_oldQ);
	if (_NRStepAccepted) return true;
	if (_Q > _oldQ) {
		_NRStepAccepted = true;
		return true;
	}

	// else
	_Q              = _oldQ;
	for (const auto f : _functions) f->rejectProposedParameters();

	return false;
}

void TModelRecal::adjustParametersPostEstimation() {
	for (const auto f : _functions) _intercept.intercept() += f->adjustParametersPostEstimation();
}

double TModelRecal::getSteepestGradient() const noexcept {
	if (_NRStepAccepted) return 0.0;
	double maxF = 0.0;
	for (unsigned int i = 0; i < numParameters(); ++i) maxF = std::max(maxF, fabs(_F(i)));
	return maxF;
}

} // namespace SequencingError
} // namespace GenotypeLikelihoods
