/*
 * TRecalibrationEMModule.cpp
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#include "TSequencingErrorCovariateFunction.h"

#include <algorithm>
#include <iterator>

#include "commonWeakTypes.h"
#include "mathFunctions.h"
#include "probability.h"
#include "stringFunctions.h"
#include "toString.h"
#include "TProbabilityDistributions.h"

#include "devtools.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

using coretools::str::toString;

namespace impl {
double normalizeParameters(std::vector<double> &betas) noexcept {
	const double mean = std::accumulate(betas.begin(), betas.end(), 0.) / betas.size();
	for (auto &bi : betas) { bi -= mean; }
	return mean;
}
} // namespace impl

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction
//--------------------------------------------------------------

void TFunction::_initializeValues(const std::vector<std::string> &betas) {
	if (!betas.empty()) {
		if (betas.size() != numParameters()) {
			throw toString("Failed to initialize recalibration module: wrong number of values (", betas.size(),
						   " instead of ", numParameters(), ")!");
		}
		std::transform(betas.cbegin(), betas.cend(), _begin(), coretools::str::convertStringCheck<double>);
	}
}

void TFunction::proposeNewParameters(const arma::mat &JxF, uint16_t &index, double lambda) noexcept {
	// update new ones
	std::copy(_cbegin(), _cend(), _obegin());
	for (auto it = _begin(); it != _end(); ++it) {
		*it -= lambda * JxF(index++);
	}
}

void TFunction::rejectProposedParameters() noexcept {
	std::copy(_obegin(), _oend(), _begin());
}

std::string TFunction::modelString() const {
	std::string s = typeString() + "[" + std::accumulate(_cbegin(), _cend(), std::string{}, [](auto tot, auto b){return tot + toString(b) + ",";});
	return s.substr(0, s.size() - 1) + "]";
}

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_probit
//--------------------------------------------------------------
TProbit::TProbitTmpStorage::TProbitTmpStorage(const std::array<double, 3> &betas, uint16_t q) {
	using namespace coretools::probdist;
	const double z        = betas[1] + betas[2] * q;
	cumulDens_Phi         = coretools::probdist::TNormalDistr::cumulativeDensity(z, 0, 1);
	normalDens_phi        = coretools::probdist::TNormalDistr::density(z, 0, 1);
	eta                   = cumulDens_Phi * betas[0];
	normalDens_q          = normalDens_phi * q;
	normalDens_Beta1      = normalDens_phi * betas[0];
	normalDens_Beta1_q    = normalDens_Beta1 * q;
	normalDens_Beta1_z    = normalDens_Beta1 * z;
	normalDens_Beta1_q_z  = normalDens_Beta1_q * z;
	normalDens_Beta1_q2_z = normalDens_Beta1_q_z * q;
}

void TProbit::_expandTmpStorage(uint16_t MaxValue) const {
	for (uint16_t q = _tmpStorage.size(); q <= MaxValue; ++q) { _tmpStorage.emplace_back(_betas, q); }
}

TProbit::TProbit(uint16_t FirstParameterIndex, const std::vector<std::string> &betas)
	: TFunction(FirstParameterIndex) {
	_initializeValues(betas);
	_expandTmpStorage(128);
}

double TProbit::getEtaTerm(uint16_t val) const noexcept {
	if (val >= _tmpStorage.size()) { _expandTmpStorage(val); }
	return _tmpStorage[val].eta;
}

T1stDerivative TProbit::get1stDerivatives(uint16_t val, size_t i) const noexcept {
	if (val >= _tmpStorage.size()) { _expandTmpStorage(val); }
	switch (i) {
	case 0: return {firstParameterIndex(), _tmpStorage[val].cumulDens_Phi};
	case 1: return {firstParameterIndex() + 1, _tmpStorage[val].normalDens_Beta1};
	default: return {firstParameterIndex() + 1, _tmpStorage[val].normalDens_Beta1_q};
	}
}

T2ndDerivative TProbit::get2ndDerivatives(uint16_t val, size_t i, size_t j) const noexcept {
	if (val >= _tmpStorage.size()) { _expandTmpStorage(val); }
	switch (i) {
	case 0:
		switch (j) {
		case 0: return {firstParameterIndex(), firstParameterIndex(), 0.};
		case 1: return {firstParameterIndex(), firstParameterIndex() + 1, -_tmpStorage[val].normalDens_phi};
		default: return {firstParameterIndex(), firstParameterIndex() + 2, -_tmpStorage[val].normalDens_q};
		}
	case 1:
		switch (j) {
		case 0: return {firstParameterIndex() + 1, firstParameterIndex(), 0.};
		case 1: return {firstParameterIndex() + 1, firstParameterIndex() + 1, -_tmpStorage[val].normalDens_Beta1_z};
		default: return {firstParameterIndex() + 1, firstParameterIndex() + 2, -_tmpStorage[val].normalDens_Beta1_q_z};
		}
	default: // 2
		switch (j) {
		case 0: return {firstParameterIndex() + 2, firstParameterIndex(), 0.};
		case 1: return {firstParameterIndex() + 2, firstParameterIndex() + 1, 0.};
		default: return {firstParameterIndex() + 2, firstParameterIndex() + 2, -_tmpStorage[val].normalDens_Beta1_q2_z};
		}
	}
}

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specific
//--------------------------------------------------------------
TSpecific::TSpecific(uint16_t FirstParameterIndex, const std::vector<std::string> &betas)
	: TFunction(FirstParameterIndex) {
	_resize(betas.size());
	_initializeValues(betas);
}

void TSpecific::_resize(uint16_t size) {
	_betas.resize(size);
	_oldBetas.resize(size);
}

void TSpecific::_adjustValueRanges(const std::vector<uint16_t> &values) {
	// initialize with maximum
	_resize(*std::max_element(values.begin(), values.end()) + 1);

	// check that each value from 0 to max is actually used!
	std::vector<bool> found(numParameters(), false);
	for (auto &i : values) { found[i] = true; }
	if (const auto f = std::find(found.cbegin(), found.cend(), false); f != found.cend()) {
		throw "Can not adjust value range for recal function '" + name + "': value " +
			toString(std::distance(f, found.cbegin())) + " is < max value but never used." +
			"\nConsider using recal function '" + TSpecificMap::name + "'.";
	}
}

double TSpecific::adjustParametersPostEstimation() noexcept { return impl::normalizeParameters(_betas); }

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specificMap
//--------------------------------------------------------------
void TSpecificMap::_resize(size_t NumParameters) {
	_betas.resize(NumParameters);
	_oldBetas.resize(NumParameters);
}

void TSpecificMap::_initMapFromVector(const std::vector<uint16_t> &values) {
	// find largest value
	const auto max = std::max((uint16_t)_indexMap.size(), *std::max_element(values.begin(), values.end()));

	// create map
	_indexMap.clear();
	_indexMap.resize(max + 1);

	for (size_t i = 0; i < values.size(); ++i) {
		_indexMap[values[i]].index = i;
		_indexMap[values[i]].index = true;
	}
}

TSpecificMap::TSpecificMap(uint16_t FirstParameterIndex, const std::vector<std::string> &betas)
	: TFunction(FirstParameterIndex) {
	// parse values as pairs separated by a colon (:)
	std::vector<uint16_t> values;
	for (std::string s : betas) {
		size_t pos = s.find(':');
		if (pos == std::string::npos) { throw "Can not parse value '" + s + "': missing ':'!"; }
		uint16_t val = coretools::str::convertStringCheck<uint16_t>(s.substr(0, pos));
		if (std::find(values.begin(), values.end(), val) != values.end()) {
			throw "Duplicate entry for key " + toString(val) + "!";
		}
		values.push_back(val);
		_betas.push_back(coretools::str::convertStringCheck<double>(s.substr(pos + 1)));
	}
	// init map
	_initMapFromVector(values);
}

void TSpecificMap::_adjustValueRanges(const std::vector<uint16_t> &values) {
	_resize(values.size());
	_initMapFromVector(values);
}

T1stDerivative TSpecificMap::get1stDerivatives(uint16_t val, size_t) const noexcept {
	return {firstParameterIndex() + _indexMap[val].index, 1.0};
}

double TSpecificMap::adjustParametersPostEstimation() noexcept { return impl::normalizeParameters(_betas); }

} // namespace SequencingError
} // namespace GenotypeLikelihoods
