/*
 * TRecalibrationEMModule.cpp
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */


#include "../GenotypeLikelihoods/TSequencingErrorCovariateFunction.h"

namespace GenotypeLikelihoods{

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction
//--------------------------------------------------------------
void TSequencingErrorCovariateFunction::_init(const uint16_t FirstParameterIndex){
	_moduleName = SequencingErrorCovariateFunction_none;
	_numParameters = 0;
	_firstParameterIndex = FirstParameterIndex;
	_numNonZeroFirstDerivatives = 0;
	_numNonZeroSecondDerivatives = 0;
	doTransformation = false;
	transformationMap = nullptr;
};

void TSequencingErrorCovariateFunction::_initializeBetas(){
	_betas.resize(_numParameters);
	_oldBetas.resize(_numParameters);

	for(uint16_t i = 0; i < _numParameters; ++i){
		_betas[i] = 0.0;
		_oldBetas[i] = 0.0;
	}
};

void TSequencingErrorCovariateFunction::_initializValues(std::vector<std::string> & values){
	if(!values.empty()){
		if(values.size() != _numParameters){
			throw "Failed to initialize recalibration module: wrong number of values (" + toString(values.size()) + " instead of " + toString(_numParameters) + ")!";
		}

		for(size_t i=0; i<values.size(); ++i){
			_betas[i] = stringToDoubleCheck(values[i]);
		}
	}
};

double TSequencingErrorCovariateFunction::_getAsDouble(const uint16_t val) const{
	if(doTransformation){
		return transformationMap->get(val);
	} else {
		return (double) val;
	}
};

double TSequencingErrorCovariateFunction::_normalizeParameters(){
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

bool TSequencingErrorCovariateFunction::checkValueRange(std::vector<uint16_t> & values){
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


void TSequencingErrorCovariateFunction::proposeNewParameters(const arma::mat & JxF, uint16_t & index, double & lambda){
	//save old parameters
	for(uint16_t i=0; i<_numParameters; ++i)
		_oldBetas[i] = _betas[i];

	//update new ones
	for(uint16_t i=0; i<_numParameters; ++i){
		_betas[i] = _oldBetas[i] - lambda * JxF(index);
		++index;
	}
};

std::string TSequencingErrorCovariateFunction::getModelString() const{
	std::string s = _moduleName + "[";
	for(size_t i=0; i<_numParameters; ++i){
		if(i>0){
			s += ",";
		}
		s += toString(_betas[i]);
	}
	s += "]";
	return s;
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_intercept
//--------------------------------------------------------------
void TSequencingErrorCovariateFunction_intercept::_init(){
	_moduleName = SequencingErrorCovariateFunction_intercept;
	_numParameters = 1;
	_numNonZeroFirstDerivatives = 1;
	_numNonZeroSecondDerivatives = 0;
	_initializeBetas();
};

TSequencingErrorCovariateFunction_intercept::TSequencingErrorCovariateFunction_intercept(){
	_init();
}

TSequencingErrorCovariateFunction_intercept::TSequencingErrorCovariateFunction_intercept(const uint16_t FirstParameterIndex):TSequencingErrorCovariateFunction(FirstParameterIndex){
	_init();
};

TSequencingErrorCovariateFunction_intercept::TSequencingErrorCovariateFunction_intercept(const uint16_t FirstParameterIndex, std::vector<std::string> & values):TSequencingErrorCovariateFunction(FirstParameterIndex){
	_init();
	_initializValues(values);
};

void TSequencingErrorCovariateFunction_intercept::initialize(const uint16_t FirstParameterIndex, std::vector<std::string> & values){
	_initializValues(values);
};

void TSequencingErrorCovariateFunction_intercept::setIntercept(const double val){
	_betas[0] = val;
};

void TSequencingErrorCovariateFunction_intercept::addToIntercept(const double val){
	_betas[0] += val;
};

double TSequencingErrorCovariateFunction_intercept::getIntercept() const{
	return _betas[0];
};

double TSequencingErrorCovariateFunction_intercept::getEtaTerm() const{
	return _betas[0];
};

double TSequencingErrorCovariateFunction_intercept::getEtaTerm(const uint16_t val){
	return _betas[0];
};

void TSequencingErrorCovariateFunction_intercept::fillDerivatives(TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const{
	first.add(_firstParameterIndex, 1.0);
};

void TSequencingErrorCovariateFunction_intercept::fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second){
	first.add(_firstParameterIndex, 1.0);
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_polynomial
//--------------------------------------------------------------
void TSequencingErrorCovariateFunction_polynomial::_init(const size_t order){
	if(order < 1)
		throw "Order of polynomial covariate function must be at least 1!";
	_moduleName = SequencingErrorCovariateFunction_polynomial;
	_numParameters = order;
	_numNonZeroFirstDerivatives = order;
	_numNonZeroSecondDerivatives = 0;
	_initializeBetas();
};

TSequencingErrorCovariateFunction_polynomial::TSequencingErrorCovariateFunction_polynomial(const uint16_t FirstParameterIndex, const size_t order):TSequencingErrorCovariateFunction(FirstParameterIndex){
	_init(order);
};

TSequencingErrorCovariateFunction_polynomial::TSequencingErrorCovariateFunction_polynomial(const uint16_t FirstParameterIndex, std::vector<std::string> & values):TSequencingErrorCovariateFunction(FirstParameterIndex){
	_init(values.size());
	_initializValues(values);
};

double TSequencingErrorCovariateFunction_polynomial::getEtaTerm(const uint16_t val) const{
	double valAsDouble = _getAsDouble(val);
	double tmp = valAsDouble;
	double sum = _betas[0] * tmp;
	for(size_t i=1; i<_numParameters; ++i){
		tmp *= valAsDouble;
		sum += _betas[i] * tmp;
	}

	return sum;
};

void TSequencingErrorCovariateFunction_polynomial::fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const{
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
TProbitTmpStorage::TProbitTmpStorage(const std::vector<double> & betas, const uint16_t q){
	double z = betas[1] + betas[2] * (double) q;
	_cumulDens_Phi[q] = TNormalDistr::cumulativeDistrFunction(z);
	_normalDens_phi[q] = TNormalDistr::density(z);
	_eta[q] = _cumulDens_Phi[q] * betas[0];
	_normalDens_q[q] = _normalDens_phi[q] * (double) q;
	_normalDens_Beta1[q] = _normalDens_phi[q] * betas[0];
	_normalDens_Beta1_q[q] = _normalDens_Beta1[q] * (double) q;
	_normalDens_Beta1_z[q] = _normalDens_Beta1[q] * z;
	_normalDens_Beta1_q_z[q] = _normalDens_Beta1_q[q] * z;
	_normalDens_Beta1_q2_z[q] = _normalDens_Beta1_q_z[q] * (double) q;
};

void TRecalibrationEMCovariateFunction_probit::_init(const uint16_t MaxValue){
	_moduleName = SequencingErrorCovariateFunction_probit;
	_numParameters = 3;
	_numNonZeroFirstDerivatives = 3;
	_numNonZeroSecondDerivatives = 6;
	_initializeBetas();
	_secondParameterIndex = _firstParameterIndex + 1;
	_thirdParameterIndex = _firstParameterIndex + 2;

	//prepare tmp storage for phi and Phi
	if(MaxValue < 1){
		MaxValue = 128;
	}
	_expandTmpStorage(MaxValue);
};

void TRecalibrationEMCovariateFunction_probit::_expandTmpStorage(const uint16_t MaxValue) const{
	for(uint16_t q = _maxValue + 1; q <= MaxValue; ++q){
		_tmpStorage.emplace_back(_betas, q);
	}
	_maxValue = MaxValue;
};

TRecalibrationEMCovariateFunction_probit::TRecalibrationEMCovariateFunction_probit(const uint16_t FirstParameterIndex, const uint16_t MaxValue):TSequencingErrorCovariateFunction(FirstParameterIndex){
	_init(MaxValue);
};

TRecalibrationEMCovariateFunction_probit::TRecalibrationEMCovariateFunction_probit(const uint16_t FirstParameterIndex, std::vector<std::string> & values):TSequencingErrorCovariateFunction(FirstParameterIndex){
	_init(0);
	_initializValues(values);
};

double TRecalibrationEMCovariateFunction_probit::getEtaTerm(const uint16_t val) const{
	if(val > _maxValue){
		_expandTmpStorage(val);
	}

	return _tmpStorage[val]._eta;
};

void TRecalibrationEMCovariateFunction_probit::fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const{
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
TSequencingErrorCovariateFunction_specific::TSequencingErrorCovariateFunction_specific(const uint16_t FirstParameterIndex, const uint16_t MaxValue):TSequencingErrorCovariateFunction(FirstParameterIndex){
	_init(MaxValue);
};

TSequencingErrorCovariateFunction_specific::TSequencingErrorCovariateFunction_specific(const uint16_t FirstParameterIndex, std::vector<std::string> & values):TSequencingErrorCovariateFunction(FirstParameterIndex){
	//are there any values provided?
	if(values.size() < 1){
		throw "Failed to initialize TRecalibrationEMCovariateFunction_specificMap: values are empty!";
	}

	//init
	_init(values.size() - 1);

	//now copy values
	_initializValues(values);
};

void TSequencingErrorCovariateFunction_specific::_init(const uint16_t MaxValue){
	_moduleName = SequencingErrorCovariateFunction_specific;
	_maxValue = MaxValue;
	_numParameters = _maxValue + 1;
	_numNonZeroFirstDerivatives = 1;
	_numNonZeroSecondDerivatives = 0;
	_initializeBetas();
};

void TSequencingErrorCovariateFunction_specific::fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const{
	first.add(_firstParameterIndex + val, 1.0);
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specificMap
//--------------------------------------------------------------
void TSequencingErrorCovariateFunction_specificMap::_init(size_t NumParameters){
	_moduleName = SequencingErrorCovariateFunction_specific;
	_numParameters = NumParameters;
	_numNonZeroFirstDerivatives = 1;
	_numNonZeroSecondDerivatives = 0;
	_initializeBetas();
};

void TSequencingErrorCovariateFunction_specificMap::_createIndexMap(){
	_mapSize = _maxValue + 1;
	_indexMap = new uint16_t[_mapSize];
	_valueUsed = new bool[_mapSize];
	for(uint16_t i=0; i<_mapSize; ++i){
		_indexMap[i] = 0;
		_valueUsed[i] = false;
	}
};

void TSequencingErrorCovariateFunction_specificMap::_freeIndexMap(){
	delete[] _indexMap;
	delete[] _valueUsed;
};

TSequencingErrorCovariateFunction_specificMap::TSequencingErrorCovariateFunction_specificMap(const uint16_t FirstParameterIndex, std::vector<uint16_t> & values):TSequencingErrorCovariateFunction(FirstParameterIndex){
	_init(values.size());

	//find largest value
	_maxValue = 0;
	for(uint16_t val : values){
		if(val > _maxValue)
			_maxValue = val;
	}

	//create map
	_createIndexMap();
	for(size_t i=0; i<values.size(); ++i){
		_indexMap[values[i]] = i;
		_valueUsed[values[i]] = true;
	}
};

TSequencingErrorCovariateFunction_specificMap::TSequencingErrorCovariateFunction_specificMap(const uint16_t FirstParameterIndex, std::vector<std::string> & values):TSequencingErrorCovariateFunction(FirstParameterIndex){
	//are there any values provided?
	if(values.size() < 1){
		throw "Failed to initialize TRecalibrationEMCovariateFunction_specificMap: values are empty!";
	}

	//init
	_init(values.size());

	//parse values as pairs separated by a colon (:)
	std::map<uint16_t, double> tmp;
	for(std::string s : values){
		size_t pos = s.find(':');
		if(pos == std::string::npos){
			throw "Can not parse value '" + s + "': missing ':'!";
		}
		uint16_t key = stringToIntCheck(s.substr(0, pos));
		if(tmp.find(key) != tmp.end()){
			throw "Duplicate entry for key " + toString(key) + "!";
		}
		tmp.emplace(key, stringToDoubleCheck(s.substr(pos+1)));
	}

	//find largest value
	_maxValue = 0;
	for(auto it = tmp.begin(); it != tmp.end(); ++it){
		if(it->first > _maxValue)
			_maxValue = it->first;
	}

	//create map
	_createIndexMap();
	for(auto it = tmp.begin(); it != tmp.end(); ++it){
		_indexMap[it->first] = it->second;
		_valueUsed[it->first] = true;
	}
};

void TSequencingErrorCovariateFunction_specificMap::fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const{
	first.add(_firstParameterIndex + _indexMap[val], 1.0);
};


}; //end namespace
