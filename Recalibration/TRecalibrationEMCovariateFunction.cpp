/*
 * TRecalibrationEMModule.cpp
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */


#include <TRecalibrationEMCovariateFunction.h>


//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction
//--------------------------------------------------------------
void TRecalibrationEMCovariateFunction::_freeBetas(){
	if(_initialized){
		delete[] _betas;
		delete[] _oldBetas;
	}
};

void TRecalibrationEMCovariateFunction::_initializeBetas(){
	_freeBetas();
	_betas = new double[_numParameters];
	_oldBetas = new double[_numParameters];
};

void TRecalibrationEMCovariateFunction::_initializValues(std::vector<std::string> & values){
	if(!values.empty()){
		if(values.size() != _numParameters){
			throw "Failed to initialize recalibration module: wrong number of values (" + toString(values.size()) + " instead of " + toString(_numParameters) + ")!";
		}

		for(size_t i=0; i<values.size(); ++i){
			_betas[i] = values[i];
		}
	}
};

double TRecalibrationEMCovariateFunction::_getAsDouble(const uint16_t val){
	if(doTransformation){
		return transformationMap->get(val);
	} else {
		return (double) val;
	}
};

bool TRecalibrationEMCovariateFunction::checkValueRange(std::vector<uint16_t> & values){
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


void TRecalibrationEMCovariateFunction::proposeNewParameters(const arma::mat & JxF, size_t & index, double & lambda){
	//save old parameters
	for(size_t i=0; i<_numParameters; ++i)
		_oldBetas[i] = _betas[i];

	//update new ones
	for(size_t i=0; i<_numParameters; ++i){
		_betas[i] = _oldBetas[i] - lambda * JxF(index);
		++index;
	}
};

std::string TRecalibrationEMCovariateFunction::getModelString(){
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
void TRecalibrationEMCovariateFunction_intercept::_init(){
	_moduleName = RecalModuleFunctionName_intercept;
	_numParameters = 1;
	_numNonZeroFirstDerivatives = 1;
	_numNonZeroSecondDerivatives = 0;
	_initializeBetas();
};

TRecalibrationEMCovariateFunction_intercept::TRecalibrationEMCovariateFunction_intercept(const size_t FirstParameterIndex):TRecalibrationEMCovariateFunction(FirstParameterIndex){
	_init();
};

TRecalibrationEMCovariateFunction_intercept::TRecalibrationEMCovariateFunction_intercept(const size_t FirstParameterIndex, std::vector<std::string> & values):TRecalibrationEMCovariateFunction(FirstParameterIndex){
	_moduleName = RecalModuleFunctionName_intercept;
	_init();
	_initializValues(values);
};

void TRecalibrationEMCovariateFunction_intercept::fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second){
	first.add(_firstParameterIndex, 1.0);
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_polynomial
//--------------------------------------------------------------
void TRecalibrationEMCovariateFunction_polynomial::_init(const size_t order){
	_moduleName = RecalModuleFunctionName_polynomial;
	_numParameters = order;
	_numNonZeroFirstDerivatives = order;
	_numNonZeroSecondDerivatives = 0;
	_initializeBetas();
};


TRecalibrationEMCovariateFunction_polynomial::TRecalibrationEMCovariateFunction_polynomial(const size_t FirstParameterIndex, int order):TRecalibrationEMCovariateFunction(FirstParameterIndex){
	_init(order);
};

TRecalibrationEMCovariateFunction_polynomial::TRecalibrationEMCovariateFunction_polynomial(const size_t FirstParameterIndex, std::vector<std::string> & values):TRecalibrationEMCovariateFunction(FirstParameterIndex){
	_init(values.size());
	_initializValues(values);
};

void TRecalibrationEMCovariateFunction_polynomial::fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second){
	double tmp = (double) val;
	first.add(_firstParameterIndex, tmp);
	for(size_t i=1; i<_numParameters; ++i){
		tmp *= (double) val;
		first.add(_firstParameterIndex + i, tmp);
	}
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specific
//--------------------------------------------------------------
TRecalibrationEMCovariateFunction_specific::TRecalibrationEMCovariateFunction_specific(const size_t FirstParameterIndex, size_t MaxValue):TRecalibrationEMCovariateFunction(FirstParameterIndex){
	_init(MaxValue);
};

TRecalibrationEMCovariateFunction_specific::TRecalibrationEMCovariateFunction_specific(const size_t FirstParameterIndex, std::vector<std::string> & values):TRecalibrationEMCovariateFunction(FirstParameterIndex){
	//are there any values provided?
	if(values.size() < 1){
		throw "Failed to initialize TRecalibrationEMCovariateFunction_specificMap: values are empty!";
	}

	//init
	_init(values.size() - 1);

	//now copy values
	_initializValues(values);
};

void TRecalibrationEMCovariateFunction_specific::_init(const size_t MaxValue){
	_moduleName = RecalModuleFunctionName_specific;
	_maxValue = MaxValue;
	_numParameters = _maxValue + 1;
	_numNonZeroFirstDerivatives = 1;
	_numNonZeroSecondDerivatives = 0;
	_initializeBetas();
};

void TRecalibrationEMCovariateFunction_specific::fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second){
	first.add(_firstParameterIndex + val, 1.0);
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specificMap
//--------------------------------------------------------------
void TRecalibrationEMCovariateFunction_specificMap::_init(size_t NumParameters){
	_moduleName = RecalModuleFunctionName_specific;
	_numParameters = NumParameters;
	_numNonZeroFirstDerivatives = 1;
	_numNonZeroSecondDerivatives = 0;
	_initializeBetas();
};

void TRecalibrationEMCovariateFunction_specificMap::_createIndexMap(){
	_mapSize = _maxValue + 1;
	_indexMap = new uint16_t[_mapSize];
	_valueUsed = new bool[_mapSize];
	for(uint16_t i=0; i<_mapSize; ++i){
		_indexMap[i] = 0;
		_valueUsed[i] = false;
	}
};

void TRecalibrationEMCovariateFunction_specificMap::_freeIndexMap(){
	delete[] _indexMap;
	delete[] _valueUsed;
};

TRecalibrationEMCovariateFunction_specificMap::TRecalibrationEMCovariateFunction_specificMap(const size_t FirstParameterIndex, std::vector<uint16_t> & values):TRecalibrationEMCovariateFunction(FirstParameterIndex){
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

TRecalibrationEMCovariateFunction_specificMap::TRecalibrationEMCovariateFunction_specificMap(const size_t FirstParameterIndex, std::vector<std::string> & values):TRecalibrationEMCovariateFunction(FirstParameterIndex){
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

void TRecalibrationEMCovariateFunction_specificMap::fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second){
	first.add(_firstParameterIndex + _indexMap[val], 1.0);
};
