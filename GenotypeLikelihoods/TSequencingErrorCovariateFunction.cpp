/*
 * TRecalibrationEMModule.cpp
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */


#include "TSequencingErrorCovariateFunction.h"
#include "stringFunctions.h"
#include <algorithm>

namespace GenotypeLikelihoods{
namespace SequencingError {

namespace /*anonymous*/ {
void initializValues(TFunction &f, const std::vector<std::string> &values) {
	if (!values.empty()) {
		if (values.size() != f.numParameters()) {
			throw coretools::str::toString("Failed to initialize recalibration module: wrong number of values (",
						       values.size(), " instead of ", f.numParameters(), ")!");
		}

		for (size_t i = 0; i < values.size(); ++i) {
			f.beta(i) = coretools::str::convertStringCheck<double>(values[i]);
		}
	}
}

double normalizeParameters(std::vector<double> &betas) noexcept {
	double mean = 0.0;
	for (const auto bi : betas) { mean += bi; }
	mean /= betas.size();
	for (auto &bi : betas) { bi -= mean; }
	return mean;
}
} // namespace

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction
//--------------------------------------------------------------

bool TFunction::checkValueRange(const std::vector<uint16_t> & values) const noexcept{
	//check range for each value
	for(uint16_t val : values){
		if(!checkValueRange(val))
			return false;
	}
	return true;
}

void TFunction::proposeNewParameters(const arma::mat & JxF, uint16_t & index, double & lambda) noexcept {
	//update new ones
	for (uint16_t i = 0; i < numParameters(); ++i) {
		_oldBeta(i) = beta(i);
		beta(i)    -= lambda * JxF(index);
		++index;
	}
}

void TFunction::rejectProposedParameters() noexcept {
	for (size_t i = 0; i < numParameters(); ++i) beta(i) = _oldBeta(i);
}

std::string TFunction::getModelString() const{
	using coretools::str::toString;
	std::string s = typeString() + "[" + toString(beta(0));
	for (size_t i = 1; i < numParameters(); ++i) s += "," + toString(beta(0));
	return s + "]";
}

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_intercept
//--------------------------------------------------------------

TIntercept::TIntercept(uint16_t FirstParameterIndex, const std::string &value)
	: TFunction(FirstParameterIndex) {
	_beta = coretools::str::convertStringCheck<double>(value);
}

void TIntercept::fillDerivatives(TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives &) const noexcept {
	first.add(firstParameterIndex(), 1.0);
}

void TIntercept::fillDerivatives(uint16_t , TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives &) const noexcept {
	first.add(firstParameterIndex(), 1.0);
}

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_polynomial
//--------------------------------------------------------------
void TPolynomial::_init(size_t order){
	if(order < 1)
		throw "Order of polynomial covariate function must be at least 1!";
	_order = order;
	_betas.resize(order);
	_oldBetas.resize(order);
}

	TPolynomial::TPolynomial(uint16_t FirstParameterIndex, size_t order, TRecalibrationEMTransformationMap* transformationMap):TFunction(FirstParameterIndex), _transformationMap(transformationMap){
	_init(order);
}

	TPolynomial::TPolynomial(uint16_t FirstParameterIndex, const std::vector<std::string> & values, TRecalibrationEMTransformationMap* transformationMap):TFunction(FirstParameterIndex), _transformationMap(transformationMap){
	_init(values.size());
	initializValues(*this, values);
}

double TPolynomial::getEtaTerm(uint16_t val) const noexcept {
	const double v = _getAsDouble(val);
	double vpi = v;
	double sum = _betas[0] * vpi;

	for(size_t i=1; i<numParameters(); ++i){
		vpi *= v;
		sum += _betas[i] * vpi;
	}
	return sum;
}

void TPolynomial::fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives &) const noexcept {
	const double v = _getAsDouble(val);
	double vpi     = v;
	first.add(firstParameterIndex(), vpi);
	for(size_t i=1; i<numParameters(); ++i){
		vpi *= v;
		first.add(firstParameterIndex() + i, vpi);
	}
}

bool TPolynomial::checkValueRange(uint16_t val) const noexcept {
	return _transformationMap ? _transformationMap->checkRange(val) : true;
}

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_probit
//--------------------------------------------------------------
TProbit::TProbitTmpStorage::TProbitTmpStorage(const std::array<double, 3> & betas, uint16_t q){
	const double z = betas[1] + betas[2] * q;
	cumulDens_Phi = coretools::TNormalDistr::cumulativeDistrFunction(z);
	normalDens_phi = coretools::TNormalDistr::density(z);
	eta = cumulDens_Phi * betas[0];
	normalDens_q = normalDens_phi * q;
	normalDens_Beta1 = normalDens_phi * betas[0];
	normalDens_Beta1_q = normalDens_Beta1 * q;
	normalDens_Beta1_z = normalDens_Beta1 * z;
	normalDens_Beta1_q_z = normalDens_Beta1_q * z;
	normalDens_Beta1_q2_z = normalDens_Beta1_q_z * q;
}

void TProbit::_init(uint16_t MaxValue){
	//prepare tmp storage for phi and Phi
	_expandTmpStorage(MaxValue < 1 ? 128 : MaxValue);
}

void TProbit::_expandTmpStorage(uint16_t MaxValue) const{
	for(uint16_t q = _maxValue + 1; q <= MaxValue; ++q){
		_tmpStorage.emplace_back(_betas, q);
	}
	_maxValue = MaxValue;
}

TProbit::TProbit(uint16_t FirstParameterIndex, uint16_t MaxValue):TFunction(FirstParameterIndex){
	_init(MaxValue);
}

TProbit::TProbit(uint16_t FirstParameterIndex, const std::vector<std::string> & values):TFunction(FirstParameterIndex){
	_init(0);
	initializValues(*this, values);
}

double TProbit::getEtaTerm(uint16_t val) const noexcept{
	if(val > _maxValue){
		_expandTmpStorage(val);
	}
	return _tmpStorage[val].eta;
}

void TProbit::fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const noexcept {
	if(val > _maxValue){
		_expandTmpStorage(val);
	}
	const auto i1 = firstParameterIndex();
	const auto i2 = i1 + 1;
	const auto i3 = i1 + 2;

	//first derivatives for each parameter
	first.add(i1, _tmpStorage[val].cumulDens_Phi);
	first.add(i2, _tmpStorage[val].normalDens_Beta1);
	first.add(i3, _tmpStorage[val].normalDens_Beta1_q);

	//add second derivatives
	second.add(i2, i2, - _tmpStorage[val].normalDens_Beta1_z);
	second.add(i3, i3, - _tmpStorage[val].normalDens_Beta1_q2_z);

	second.add(i1, i2, _tmpStorage[val].normalDens_phi);
	second.add(i1, i3, _tmpStorage[val].normalDens_q);
	second.add(i2, i3, - _tmpStorage[val].normalDens_Beta1_q_z);
}

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specific
//--------------------------------------------------------------
TSpecific::TSpecific(uint16_t FirstParameterIndex, uint16_t MaxValue):TFunction(FirstParameterIndex){
	_init(MaxValue);
}

TSpecific::TSpecific(uint16_t FirstParameterIndex, const std::vector<std::string> & betas):TFunction(FirstParameterIndex){
	//init
	_init(betas.size() - 1);

	//now copy values
	initializValues(*this, betas);
}

void TSpecific::_init(uint16_t MaxValue){
	_maxValue = MaxValue;
	_betas.reserve(numParameters());
	_oldBetas.reserve(numParameters());
}

void TSpecific::adjustValueRanges(const std::vector<uint16_t> & usedValues){
	//initialize with maximum
	uint16_t max = *std::max_element(usedValues.begin(), usedValues.end());
	_init(max);

	//check that each value from 0 to max is actually used!
	std::vector<bool> found(max+1, false);
	for(auto& i : usedValues){
		found[i] = true;
	}
	for(uint16_t i=0; i<found.size(); ++i){
		if(!found[i]){
			throw "Can not adjust value range for recal function '" + name + "': value " + coretools::str::toString(i) + " is < max value but never used."
				+ "\nConsider using recal function '" + TSpecificMap::name + "'.";
		}
	}
}

void TSpecific::fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives &) const noexcept {
	first.add(firstParameterIndex() + val, 1.0);
}

double TSpecific::adjustParametersPostEstimation() noexcept { return normalizeParameters(_betas); }

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specificMap
//--------------------------------------------------------------
void TSpecificMap::_init(size_t NumParameters){
	_numParameters = NumParameters;
	_betas.reserve(numParameters());
	_oldBetas.reserve(numParameters());
}

void TSpecificMap::_initMapFromVector(const std::vector<uint16_t> & values){
	_init(values.size());

	//find largest value
	_maxValue = std::max(_maxValue, *std::max_element(values.begin(), values.end()));

	//create map
	_indexMap.clear();
	_indexMap.resize(_maxValue + 1);

	for(size_t i=0; i<values.size(); ++i){
		_indexMap[values[i]].index = i;
		_indexMap[values[i]].index = true;
	}
}

TSpecificMap::TSpecificMap(uint16_t FirstParameterIndex, const std::vector<uint16_t> & values):TFunction(FirstParameterIndex){
	_initMapFromVector(values);
}

TSpecificMap::TSpecificMap(uint16_t FirstParameterIndex, const std::vector<std::string> & values):TFunction(FirstParameterIndex){
	//parse values as pairs separated by a colon (:)
	std::vector<uint16_t> valuesUsed;
	for(std::string s : values){
		size_t pos = s.find(':');
		if(pos == std::string::npos){
			throw "Can not parse value '" + s + "': missing ':'!";
		}
		uint16_t key = coretools::str::convertStringCheck<uint16_t>(s.substr(0, pos));
		 if (std::find(valuesUsed.begin(), valuesUsed.end(), key) != valuesUsed.end()){
			throw "Duplicate entry for key " + coretools::str::toString(key) + "!";
		}
		valuesUsed.push_back(key);
		_betas.push_back(coretools::str::convertStringCheck<double>(s.substr(pos+1)));
	}
	//init map
	_initMapFromVector(valuesUsed);
}

void TSpecificMap::adjustValueRanges(const std::vector<uint16_t> & valuesUsed){
	_initMapFromVector(valuesUsed);
}

void TSpecificMap::fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives &) const noexcept {
	first.add(firstParameterIndex() + _indexMap[val].index, 1.0);
}

double TSpecificMap::adjustParametersPostEstimation() noexcept { return normalizeParameters(_betas); }

} // namespace SequencingError
} // namespace GenotypeLikelihoods
