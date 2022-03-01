/*
 * TRecalibrationEMModel.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#include "TSequencingErrorModel.h"
#include "probability.h"
#include <memory>
#include "devtools.h"

namespace GenotypeLikelihoods {
namespace SequencingError {
using coretools::Probability;

//*********************************************************
// TRecalibrationEMModelCovariateDefinition
// class to store model definition. Used when parsing files
//*********************************************************

TCovariateDefinition::TCovariateDefinition(const std::string &modelString) {
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
				_intercept = "0.0";
			} else {
				if (pos_2 == std::string::npos)
					throw "Unable to understand model string '" + modelString + "': missing ']'";
				_intercept = s.substr(pos_1 + 1, pos_2 - pos_1 - 1);
			}
		} else {
			_covariateFunctions.push_back({s.substr(0, pos), s.substr(pos + 1)});
		}
	}
}

void TCovariateDefinition::setIntercept(double Intercept) { _intercept = coretools::str::toString(Intercept); };

void TCovariateDefinition::addCovariate(const std::string &covariate, const std::string &function) {
	_covariateFunctions.push_back({covariate, function});
}

std::string TCovariateDefinition::getModelString() const {
	std::string modelString = "intercept[" + _intercept + "]";
	for (auto &it : _covariateFunctions) { modelString += ";" + it.covariate + ":" + it.function; }
	return modelString;
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
	for (size_t a = 0; a < vec.size(); ++a) {
		std::string &row = vec[a];
		coretools::str::trimString(row, "()");
		coretools::str::fillContainerFromString(row, r, ',');
		if (r.size() != 4)
			throw "Rho matrix has " + toString(r.size()) + " instead of 4 columns for row " + toString(a + 1) + "!";

		// fill
		for (size_t b = 0; b < 4; ++b) {
			if (a == b) {
				rho[a][b] = 0.0;
			} else {
				rho[a][b] = r[b];
			}
		}
	}
}

std::string TRho::getDefinition() const noexcept {
	std::string def;
	for (int a = 0; a < 4; ++a) {
		if (a > 0) def += ";";
		for (int b = 0; b < 4; ++b) {
			if (b > 0) def += ",";
			if (a != b)
				def += coretools::str::toString(rho[a][b]);
			else
				def += "-";
		}
	}
	return def;
}

void TRho::addBaseForEstimation(genometools::Base base, const TBaseLikelihoods &EMWeights) noexcept {
	using genometools::Base;
	using genometools::index;
	for (auto b = Base::min; b < Base::max; ++b) {
		if (base != b) rho[index(b)][index(base)] += EMWeights[b].get();
	}
}

void TRho::estimate() noexcept {
	for (int a = 0; a < 4; ++a) {
		rho[a][a] = 0.0;
		double d  = 0.;
		for (const auto r : rho[a]) d += r;
		for (auto &r : rho[a]) r /= d;
	}
}

//*********************************************************
// TModelNoRecal
//*********************************************************
Probability TModelNoRecal::getErrorRate(const BAM::TSequencedBase &base) const noexcept {
	if (base == genometools::Base::N) return Probability(1.0);
	return (Probability)base.originalQuality_phredInt;
}

genometools::PhredIntProbability TModelNoRecal::getPhredInt(const BAM::TSequencedBase &base) const noexcept {
	// Todo: change to maxProb() one available.
	if (base == genometools::Base::N) return genometools::PhredIntProbability(0);
	return base.originalQuality_phredInt;
}

void TModelNoRecal::fillBaseLikelihoods(const BAM::TSequencedBase &base,
					TBaseLikelihoods &baseLikelihoods) const noexcept {
	using genometools::Base;
	if (base == Base::N) {
		baseLikelihoods.reset();
	} else {
		const auto eps = static_cast<Probability>(base.originalQuality_phredInt);
		for (auto other = Base::min; other < Base::max; ++other) {
			if (other == base.base)
				baseLikelihoods[other] = eps.complement();
			else
				baseLikelihoods[other] = (1. / 3) * eps;
		}
	}
}

void TModelNoRecal::simulate(BAM::TSequencedBase &base, coretools::TRandomGenerator &RandomGenerator) const noexcept {
	using genometools::Base;
	if (base.base == Base::N) return;

	const auto eps = static_cast<Probability>(base.originalQuality_phredInt);
	if (RandomGenerator.getRand() < (1. / 3) * eps) {
		const int i = RandomGenerator.getRand(0, 3); // 3 bases to choose from
		base.base   = Base((index(base.base) + i) % 4);
	}
}

//*********************************************************
// TModelRecal
//*********************************************************

TModelRecal::TModelRecal(const TModelDefinition &modelDef) {
	const auto &covariateMap = modelDef.covariates;

	_rho = modelDef.rho;

	// create covariates
	_intercept.initialize(0, {covariateMap.intercept()});
	_functions.push_back(&_intercept);

	_numParameters = _intercept.numParameters();
	for (auto it = covariateMap.cbegin(); it != covariateMap.cend(); ++it) {
		// create function for each covariate
		if (it->covariate == TCovariate::name) {
			continue;
		} else if (it->covariate == TCovariate_quality::name) {
			_covariates.push_back(std::make_unique<TCovariate_quality>(_numParameters, it->function));
		} else if (it->covariate == TCovariate_position::name) {
			_covariates.push_back(std::make_unique<TCovariate_position>(_numParameters, it->function));
		} else if (it->covariate == TCovariate_context::name) {
			_covariates.push_back(std::make_unique<TCovariate_context>(_numParameters, it->function));
		} else if (it->covariate == TCovariate_fragmentLength::name) {
			_covariates.push_back(std::make_unique<TCovariate_fragmentLength>(_numParameters, it->function));
		} else if (it->covariate == TCovariate_mappingQuality::name) {
			_covariates.push_back(std::make_unique<TCovariate_mappingQuality>(_numParameters, it->function));
		} else {
			throw "Unknown recalibration covariate '" + it->covariate + "' with function " + it->function + "!";
		}
		_functions.push_back(_covariates.back()->getPointerToFunction());

		// add new parameters
		_numParameters += _covariates.back()->numParameters();
	}

	// prepare Newton-Raphson variables
	setNewtonRaphsonParamsToZero();
}

TModelRecal::TModelRecal(const TModelDefinition &modelDef, const RecalEstimatorTools::TRecalDataTable &dataTable) {
	const auto &covariateMap = modelDef.covariates;

	_rho = modelDef.rho;

	// create covariates
	_intercept.initialize(0, {covariateMap.intercept()});
	_functions.push_back(&_intercept);

	_numParameters = _intercept.numParameters();
	for (auto it = covariateMap.cbegin(); it != covariateMap.cend(); ++it) {
		// create function for each covariate
		if (it->covariate == TCovariate::name) {
			continue;
		} else if (it->covariate == TCovariate_quality::name) {
			_covariates.emplace_back(new TCovariate_quality(_numParameters, it->function, dataTable));
		} else if (it->covariate == TCovariate_position::name) {
			_covariates.emplace_back(new TCovariate_position(_numParameters, it->function, dataTable));
		} else if (it->covariate == TCovariate_context::name) {
			_covariates.emplace_back(new TCovariate_context(_numParameters, it->function, dataTable));
		} else if (it->covariate == TCovariate_fragmentLength::name) {
			_covariates.emplace_back(new TCovariate_fragmentLength(_numParameters, it->function, dataTable));
		} else if (it->covariate == TCovariate_mappingQuality::name) {
			_covariates.emplace_back(new TCovariate_mappingQuality(_numParameters, it->function, dataTable));
		} else {
			throw "Unknown recalibration covariate '" + it->covariate + "' with function " + it->function + "!";
		}

		// add new parameters
		_numParameters += _covariates.back()->numParameters();
		_functions.push_back(_covariates.back()->getPointerToFunction());
	}

	// prepare Newton-Raphson variables
	setNewtonRaphsonParamsToZero();
}

TModelDefinition TModelRecal::getModelDefinition() const {
	return TModelDefinition(getCovariateDefinition(), getRhoDefinition());
}

std::string TModelRecal::getCovariateDefinition() const noexcept {
	TCovariateDefinition def;
	def.setIntercept(_intercept.getIntercept());
	for (const auto &cov : _covariates) { def.addCovariate(cov->typeString(), cov->functionString()); }
	return def.getModelString();
}

//-------------------------------------------------
// functions to calculate error rates
//-------------------------------------------------

Probability TModelRecal::_calcEpsilon(double eta) const noexcept {
	if (eta > 23.03) return Probability(0.9999999999);
	if (eta < -23.03) return Probability(0.0000000001);

	return coretools::logistic(eta);
}

Probability TModelRecal::_calcErrorRate(const BAM::TSequencedBase &base) const noexcept {
	// eta = bta[0] + SUM_i f(q[i]), where the functions are implemented as covariate function

	auto eta = _intercept.getEtaTerm();
	for (const auto &cov : _covariates) eta += cov->getEtaTerm(base);

	return _calcEpsilon(eta);
}

Probability TModelRecal::getErrorRate(const BAM::TSequencedBase &base) const noexcept {
	if (base == genometools::Base::N) return Probability::highest();
	return _calcErrorRate(base);
}

genometools::PhredIntProbability TModelRecal::getPhredInt(const BAM::TSequencedBase &base) const noexcept {
	// Todo: change to maxProb() one available.
	if (base == genometools::Base::N) return genometools::PhredIntProbability(0);
	return genometools::PhredIntProbability(_calcErrorRate(base));
}

void TModelRecal::fillBaseLikelihoods(const BAM::TSequencedBase &base, TBaseLikelihoods &baseLikelihoods) const noexcept {
	using genometools::Base;
	if (base == Base::N) {
		baseLikelihoods.reset();
	} else {
		const auto e = _calcErrorRate(base);
		for (auto b = Base::min; b < Base::max; ++b) {
			baseLikelihoods[b] = e * _rho(b, base.base);
		}
		baseLikelihoods[base.base] = e.complement();
	}
}

//-------------------------------------------------
// functions to estimate rho
//-------------------------------------------------
void TModelRecal::prepareRhoEstimationFromEMWeights() noexcept { _rho.prepareEstimationFromEMWeights(); }

void TModelRecal::addBaseForRhoEstimation(BAM::TSequencedBase &base, const TBaseLikelihoods &EMWeights) noexcept {
	_rho.addBaseForEstimation(base.base, EMWeights);
}

void TModelRecal::estimateRho() noexcept { _rho.estimate(); }

//-------------------------------------------------
// functions for estimation
//-------------------------------------------------
std::string TModelRecal::checkParameterRange(RecalEstimatorTools::TRecalDataTable &DataTable) const {
	for (auto &cov : _covariates) {
		if (!cov->checkParameterRange(DataTable)) {
			 return "Function for covariate " + cov->typeString() + " does not cover full range of data";
		}
	}
	return "";
}

void TModelRecal::_initializeDerivatives() {
	// intercept
	size_t numNonZeroFirstDeriv  = _intercept.numNonZeroFirstDerivatives();
	size_t numNonZeroSecondDeriv = _intercept.numNonZeroSecondDerivatives();

	// covariates
	for (const auto &cov : _covariates) {
		numNonZeroFirstDeriv += cov->numNonZeroFirstDerivatives();
		numNonZeroSecondDeriv += cov->numNonZeroSecondDerivatives();
	}

	_firstDerivatives.resize(numNonZeroFirstDeriv);
	_secondDerivatives.resize(numNonZeroSecondDeriv);
}

void TModelRecal::setNewtonRaphsonParamsToZero() {
	_Jacobian.resize(numParameters(), numParameters());
	_F.resize(numParameters());
	_JxF.resize(numParameters(), 1);

	_Jacobian.zeros();
	_F.zeros();

	_initializeDerivatives();

	_numSitesAdded  = 0;
	_NRconverged    = false;
	_NRStepAccepted = false;
}

void TModelRecal::setQToZero() noexcept {
	if (!_NRconverged) {
		_oldQ = _Q;
		_Q    = 0.0;
	}
}

void TModelRecal::addToQ(const BAM::TSequencedBase &base, const TBaseLikelihoods &EM_weights_bbar_given_d) {
	if (!_NRconverged) {
		// get error rate
		const auto eps = getErrorRate(base);
		// calculate sum_bbar [ Ind(bbar=d)log(1-eps) + Ind(bbar!=d)log(eps) ]
		_Q += EM_weights_bbar_given_d[base.base].get() * log(eps.complement().get()) +
		      EM_weights_bbar_given_d[base.base].complement().get() * log(eps.get());
	}
}

void TModelRecal::addToFandJacobian(const BAM::TSequencedBase &base, const TBaseLikelihoods &EM_weights_bbar_given_d) {
	// 1) Calculate epsilon
	//--------------------
	const auto eps = getErrorRate(base).get();

	// 2 ) fill derivatives
	//--------------------
	_firstDerivatives.restart();
	_secondDerivatives.restart();

	// fill derivatives of intercept
	_intercept.fillDerivatives(0.0, _firstDerivatives, _secondDerivatives);

	// fill derivatives of covariates
	for (const auto &cov : _covariates) cov->fillDerivatives(base, _firstDerivatives, _secondDerivatives);

	// 3) add derivatives to F and Jacobian
	// calculate weights
	const auto weight1 = 1.0 - eps - EM_weights_bbar_given_d[base.base].get();
	const auto weight2 = (1.0 - eps) * eps;

	OUT(weight1);
	OUT(weight2);
	OUT(eps);
	OUT(EM_weights_bbar_given_d[base.base]);

	// add first derivatives
	for (auto it = _firstDerivatives.begin(); it != _firstDerivatives.end(); ++it) {
		// add to F
		_F(it->index) += weight1 * it->derivative;

		// add to J
		for (auto it2 = it; it2 != _firstDerivatives.end(); ++it2) {
			_Jacobian(it->index, it2->index) += weight2 * it->derivative * it2->derivative;
		}
	}

	// add second derivatives to Jacobian (happens to have the same weigth as F!)
	for (auto &it : _secondDerivatives) _Jacobian(it.index1, it.index2) += weight1 * it.derivative;

	++_numSitesAdded;
}

bool TModelRecal::solveJxF() {
	if (_NRconverged) return true;

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

void TModelRecal::proposeNewParameters(double &lambda) {
	if (!_NRStepAccepted) {
		uint16_t index = 0;
		for (const auto it : _functions) { it->proposeNewParameters(_JxF, index, lambda); }
	}
}

bool TModelRecal::acceptProposedParametersBasedOnQ() {
	if (_NRStepAccepted) return true;
	if (_Q > _oldQ) {
		_NRStepAccepted = true;
		return true;
	}
	_NRStepAccepted = false;
	_Q              = _oldQ;
	for (const auto it : _functions) it->rejectProposedParameters();

	return false;
}

void TModelRecal::adjustParametersPostEstimation() {
	for (const auto it : _functions) _intercept.addToIntercept(it->adjustParametersPostEstimation());
}

double TModelRecal::getSteepestGradient() const noexcept {
	if (_NRStepAccepted) return 0.0;
	double maxF = 0.0;
	for (unsigned int i = 0; i < numParameters(); ++i) maxF = std::max(maxF, fabs(_F(i)));
	return maxF;
}

} // namespace SequencingError
} // namespace GenotypeLikelihoods
