/*
 * TRecalibrationEMModule.cpp
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */


#include "TSequencingErrorCovariateFunction.h"

namespace GenotypeLikelihoods{
namespace SequencingError {

//define module names
extern const std::string SequencingErrorCovariateFunction_none = "none";
extern const std::string SequencingErrorCovariateFunction_intercept = "intercept";
extern const std::string SequencingErrorCovariateFunction_polynomial = "polynomial";
extern const std::string SequencingErrorCovariateFunction_probit = "probit";
extern const std::string SequencingErrorCovariateFunction_specific = "specific";
extern const std::string SequencingErrorCovariateFunction_map = "map";

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction
//--------------------------------------------------------------
void TCovariateFunction::_init(uint16_t FirstParameterIndex){
	_moduleName = SequencingErrorCovariateFunction_none;
	_numParameters = 0;
	_firstParameterIndex = FirstParameterIndex;
	_numNonZeroFirstDerivatives = 0;
	_numNonZeroSecondDerivatives = 0;
	doTransformation = false;
	transformationMap = nullptr;
};

void TCovariateFunction::_initializeBetas(){
	_betas.resize(_numParameters);
	_oldBetas.resize(_numParameters);

	for(uint16_t i = 0; i < _numParameters; ++i){
		_betas[i] = 0.0;
		_oldBetas[i] = 0.0;
	}
};

void TCovariateFunction::_initializValues(const std::vector<std::string> & values){
	if(!values.empty()){
		if(values.size() != _numParameters){
			throw coretools::str::toString("Failed to initialize recalibration module: wrong number of values (", values.size(), " instead of ", _numParameters, ")!");
		}

		for(size_t i=0; i<values.size(); ++i){
			_betas[i] = coretools::str::convertStringCheck<double>(values[i]);
		}
	}
};

double TCovariateFunction::_getAsDouble(uint16_t val) const{
	if(doTransformation){
		return transformationMap->get(val);
	} else {
		return (double) val;
	}
};

double TCovariateFunction::_normalizeParameters(){
	double mean = 0.0;
	for(uint16_t i=0; i<_numParameters; ++i){
		mean += _betas[i];
	}
	mean /= (double) _numParameters;


	for(uint16_t i=0; i<_numParameters; ++i){
		_betas[i] -= mean;
	}

	return mean;
};

bool TCovariateFunction::checkValueRange(const std::vector<uint16_t> & values) const{
	//check range for each value
	for(uint16_t val : values){
		if(!checkValueRange(val))
			return false;
	}
	//check range regarding transformation
	if(doTransformation){
		for(uint16_t val : values){
			if(!transformationMap->checkRange(val))
				return false;
		}
	}
	return true;
};


void TCovariateFunction::proposeNewParameters(const arma::mat & JxF, uint16_t & index, double & lambda){
	//save old parameters
	for(uint16_t i=0; i<_numParameters; ++i)
		_oldBetas[i] = _betas[i];

	//update new ones
	for(uint16_t i=0; i<_numParameters; ++i){
		_betas[i] = _oldBetas[i] - lambda * JxF(index);
		++index;
	}
};

std::string TCovariateFunction::getModelString() const{
	return _moduleName + "[" + coretools::str::concatenateString(_betas, ",") + "]";
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_intercept
//--------------------------------------------------------------
void TCovariateFunction_intercept::_init(){
	_moduleName = SequencingErrorCovariateFunction_intercept;
	_numParameters = 1;
	_numNonZeroFirstDerivatives = 1;
	_numNonZeroSecondDerivatives = 0;
	_initializeBetas();
};

TCovariateFunction_intercept::TCovariateFunction_intercept(){
	_init();
}

TCovariateFunction_intercept::TCovariateFunction_intercept(uint16_t FirstParameterIndex):TCovariateFunction(FirstParameterIndex){
	_init();
};

TCovariateFunction_intercept::TCovariateFunction_intercept(uint16_t FirstParameterIndex, const std::vector<std::string> & values):TCovariateFunction(FirstParameterIndex){
	_init();
	_initializValues(values);
};

void TCovariateFunction_intercept::initialize(uint16_t ){
	_init();
};

void TCovariateFunction_intercept::initialize(uint16_t , const std::vector<std::string> & values){
	_initializValues(values);
};

void TCovariateFunction_intercept::setIntercept(double val){
	_betas[0] = val;
};

void TCovariateFunction_intercept::addToIntercept(double val){
	_betas[0] += val;
};

double TCovariateFunction_intercept::getIntercept() const{
	return _betas[0];
};

double TCovariateFunction_intercept::getEtaTerm() const{
	return _betas[0];
};

double TCovariateFunction_intercept::getEtaTerm(uint16_t ) const{
	return _betas[0];
};

void TCovariateFunction_intercept::fillDerivatives(TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives &) const{
	first.add(_firstParameterIndex, 1.0);
};

void TCovariateFunction_intercept::fillDerivatives(uint16_t , TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives &) const{
	first.add(_firstParameterIndex, 1.0);
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_polynomial
//--------------------------------------------------------------
void TCovariateFunction_polynomial::_init(size_t order){
	if(order < 1)
		throw "Order of polynomial covariate function must be at least 1!";
	_moduleName = SequencingErrorCovariateFunction_polynomial;
	_numParameters = order;
	_numNonZeroFirstDerivatives = order;
	_numNonZeroSecondDerivatives = 0;
	_initializeBetas();
};

TCovariateFunction_polynomial::TCovariateFunction_polynomial(uint16_t FirstParameterIndex, size_t order):TCovariateFunction(FirstParameterIndex){
	_init(order);
};

TCovariateFunction_polynomial::TCovariateFunction_polynomial(uint16_t FirstParameterIndex, const std::vector<std::string> & values):TCovariateFunction(FirstParameterIndex){
	_init(values.size());
	_initializValues(values);
};

double TCovariateFunction_polynomial::getEtaTerm(uint16_t val) const{
	double valAsDouble = _getAsDouble(val);
	double tmp = valAsDouble;
	double sum = _betas[0] * tmp;

	for(size_t i=1; i<_numParameters; ++i){
		tmp *= valAsDouble;
		sum += _betas[i] * tmp;
	}

	return sum;
};

void TCovariateFunction_polynomial::fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives &) const{
	double valAsDouble = _getAsDouble(val);
	double tmp = valAsDouble;
	first.add(_firstParameterIndex, tmp);
	for(size_t i=1; i<_numParameters; ++i){
		tmp *= valAsDouble;
		first.add(_firstParameterIndex + i, tmp);
	}
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_probit
//--------------------------------------------------------------
TProbitTmpStorage::TProbitTmpStorage(const std::vector<double> & betas, uint16_t q){
	double z = betas[1] + betas[2] * (double) q;
	_cumulDens_Phi = coretools::TNormalDistr::cumulativeDistrFunction(z);
	_normalDens_phi = coretools::TNormalDistr::density(z);
	_eta = _cumulDens_Phi * betas[0];
	_normalDens_q = _normalDens_phi * (double) q;
	_normalDens_Beta1 = _normalDens_phi * betas[0];
	_normalDens_Beta1_q = _normalDens_Beta1 * (double) q;
	_normalDens_Beta1_z = _normalDens_Beta1 * z;
	_normalDens_Beta1_q_z = _normalDens_Beta1_q * z;
	_normalDens_Beta1_q2_z = _normalDens_Beta1_q_z * (double) q;
};

void TRecalibrationEMCovariateFunction_probit::_init(uint16_t MaxValue){
	_moduleName = SequencingErrorCovariateFunction_probit;
	_numParameters = 3;
	_numNonZeroFirstDerivatives = 3;
	_numNonZeroSecondDerivatives = 6;
	_initializeBetas();
	_secondParameterIndex = _firstParameterIndex + 1;
	_thirdParameterIndex = _firstParameterIndex + 2;

	//prepare tmp storage for phi and Phi
	if(MaxValue < 1){
		_maxValue = 128;
	}
	_expandTmpStorage(MaxValue);
};

void TRecalibrationEMCovariateFunction_probit::_expandTmpStorage(uint16_t MaxValue) const{
	for(uint16_t q = _maxValue + 1; q <= MaxValue; ++q){
		_tmpStorage.emplace_back(_betas, q);
	}
	_maxValue = MaxValue;
};

TRecalibrationEMCovariateFunction_probit::TRecalibrationEMCovariateFunction_probit(uint16_t FirstParameterIndex, uint16_t MaxValue):TCovariateFunction(FirstParameterIndex){
	_init(MaxValue);
};

TRecalibrationEMCovariateFunction_probit::TRecalibrationEMCovariateFunction_probit(uint16_t FirstParameterIndex, const std::vector<std::string> & values):TCovariateFunction(FirstParameterIndex){
	_init(0);
	_initializValues(values);
};

double TRecalibrationEMCovariateFunction_probit::getEtaTerm(uint16_t val) const{
	if(val > _maxValue){
		_expandTmpStorage(val);
	}

	return _tmpStorage[val]._eta;
};

void TRecalibrationEMCovariateFunction_probit::fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const{
	if(val > _maxValue){
		_expandTmpStorage(val);
	}

	//first derivatives for each parameter
	first.add(_firstParameterIndex, _tmpStorage[val]._cumulDens_Phi);
	first.add(_secondParameterIndex, _tmpStorage[val]._normalDens_Beta1);
	first.add(_thirdParameterIndex, _tmpStorage[val]._normalDens_Beta1_q);

	//add second derivatives
	second.add(_secondParameterIndex, _secondParameterIndex, - _tmpStorage[val]._normalDens_Beta1_z);
	second.add(_thirdParameterIndex, _thirdParameterIndex, - _tmpStorage[val]._normalDens_Beta1_q2_z);

	second.add(_firstParameterIndex, _secondParameterIndex, _tmpStorage[val]._normalDens_phi);
	second.add(_firstParameterIndex, _thirdParameterIndex, _tmpStorage[val]._normalDens_q);
	second.add(_secondParameterIndex, _thirdParameterIndex, - _tmpStorage[val]._normalDens_Beta1_q_z);
};


//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specific
//--------------------------------------------------------------
TCovariateFunction_specific::TCovariateFunction_specific(uint16_t FirstParameterIndex, uint16_t MaxValue):TCovariateFunction(FirstParameterIndex){
	_init(MaxValue);
};

TCovariateFunction_specific::TCovariateFunction_specific(uint16_t FirstParameterIndex, const std::vector<std::string> & betas):TCovariateFunction(FirstParameterIndex){
	//init
	_init(betas.size() - 1);

	//now copy values
	_initializValues(betas);
};

void TCovariateFunction_specific::_init(uint16_t MaxValue){
	_moduleName = SequencingErrorCovariateFunction_specific;
	_maxValue = MaxValue;
	_numParameters = _maxValue + 1;
	_numNonZeroFirstDerivatives = 1;
	_numNonZeroSecondDerivatives = 0;
	_initializeBetas();
};

bool TCovariateFunction_specific::checkValueRange(uint16_t val) const {
	if(val > _maxValue)
		return false;
	return true;
};

void TCovariateFunction_specific::adjustValueRanges(const std::vector<uint16_t> & usedValues){
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
			throw "Can not adjust value range for recal function '" + SequencingErrorCovariateFunction_specific + "': value " + coretools::str::toString(i) + " is < max value but never used."
				+ "\nConsider using recal function '" + SequencingErrorCovariateFunction_map + "'.";
		}
	}
};

void TCovariateFunction_specific::fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives &) const{
	first.add(_firstParameterIndex + val, 1.0);
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specificMap
//--------------------------------------------------------------
void TCovariateFunction_specificMap::_init(size_t NumParameters){
	_moduleName = SequencingErrorCovariateFunction_specific;
	_numParameters = NumParameters;
	_numNonZeroFirstDerivatives = 1;
	_numNonZeroSecondDerivatives = 0;
	_initializeBetas();
};

void TCovariateFunction_specificMap::_initMapFromVector(const std::vector<uint16_t> & values){
	_init(values.size());

	//find largest value
	_maxValue = 0;
	for(uint16_t val : values){
		if(val > _maxValue)
			_maxValue = val;
	}

	//create map
	_indexMap.clear();
	_indexMap.resize(_maxValue + 1);

	for(size_t i=0; i<values.size(); ++i){
		_indexMap[values[i]].index = i;
		_indexMap[values[i]].index = true;
	}
};

TCovariateFunction_specificMap::TCovariateFunction_specificMap(uint16_t FirstParameterIndex, const std::vector<uint16_t> & values):TCovariateFunction(FirstParameterIndex){
	_initMapFromVector(values);
};

TCovariateFunction_specificMap::TCovariateFunction_specificMap(uint16_t FirstParameterIndex, const std::vector<std::string> & values):TCovariateFunction(FirstParameterIndex){
	//parse values as pairs separated by a colon (:)
	std::vector<uint16_t> valuesUsed;
	std::vector<double> betas;
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
		betas.push_back(coretools::str::convertStringCheck<double>(s.substr(pos+1)));
	}

	//init map
	_initMapFromVector(valuesUsed);

	//copy betas
	_betas = betas;
};

bool TCovariateFunction_specificMap::checkValueRange(uint16_t val) const{
	if(val > _maxValue || !_indexMap[val].used){
		return false;
	}
	return true;
};

void TCovariateFunction_specificMap::adjustValueRanges(const std::vector<uint16_t> & valuesUsed){
	_initMapFromVector(valuesUsed);
};

void TCovariateFunction_specificMap::fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives &) const{
	first.add(_firstParameterIndex + _indexMap[val].index, 1.0);
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods
