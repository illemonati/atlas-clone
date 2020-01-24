/*
 * TRecalibrationEMModule.h
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#ifndef RECALIBRATION_TRECALIBRATIONEMCOVARIATEFUNCTION_H_
#define RECALIBRATION_TRECALIBRATIONEMCOVARIATEFUNCTION_H_

#include "../stringFunctions.h"
#include "TRecalibrationEMAuxiliaryTools.h"
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>


//define module names
#define RecalModuleFunctionName_none "none"
#define RecalModuleFunctionName_intercept "intercept"
#define RecalModuleFunctionName_polynomial "polynomial"
#define RecalModuleFunctionName_specific "specific"


//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction
// Base class for recal covariate functions
//--------------------------------------------------------------

class TRecalibrationEMCovariateFunction{
protected:
	size_t _numParameters;
	std::string _moduleName;
	double* _betas; //betas of the model
	double* _oldBetas; //use during estimation
	bool _initialized;

	//transform values?
	TRecalibrationEMTransformationMap* transformationMap;
	bool doTransformation;

	void _freeBetas();
	void _initializeBetas();
	void _initializValues(std::vector<std::string> & values);

	double _getAsDouble(const uint16_t val){
		if(doTransformation){
			return transformationMap->get(val);
		} else {
			return (double) val;
		}
	};

	friend class TRecalibrationEMModelNEW;

public:
	TRecalibrationEMCovariateFunction(){
		_moduleName = RecalModuleFunctionName_none;
		_initialized = false;
		_numParameters = 0;
		_betas = nullptr;
		_oldBetas = nullptr;

		doTransformation = false;
		transformationMap = nullptr;
	};

	virtual ~TRecalibrationEMCovariateFunction(){
		_freeBetas();
	};

	virtual bool checkValueRange(const uint16_t val){
		return true;
	};

	bool checkValueRange(std::vector<uint16_t> & values){
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

	int numParameters(){
		return numParameters;
	};

	void addTransformation(TRecalibrationEMTransformationMap* pointerToTransformationMap){
		doTransformation = true;
		transformationMap = pointerToTransformationMap;
	};

	virtual double getEtaTerm(const uint16_t val){ return 0.0; }

	void proposeNewParameters(const arma::mat & JxF, size_t & index, double & lambda){
		//save old parameters
		for(size_t i=0; i<_numParameters; ++i)
			_oldBetas[i] = _betas[i];

		//update new ones
		for(size_t i=0; i<_numParameters; ++i){
			_betas[i] = _oldBetas[i] - lambda * JxF(index);
			++index;
		}
	};

	void rejectProposedParameters(){
		for(unsigned int i=0; i<_numParameters; ++i)
			_betas[i] = _oldBetas[i];
	};

	std::string getModelString(){
		std::string s = _moduleName + "(";
		for(size_t i=0; i<_numParameters; ++i){
			if(i>0){
				s += ",";
			}
			s += toString(_betas[i]);
		}
		s += ")";
		return s;
	};
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_intercept
// An intercept term
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction_intercept:public TRecalibrationEMCovariateFunction{
protected:

public:
	TRecalibrationEMCovariateFunction_intercept(){
		_moduleName = RecalModuleFunctionName_intercept;
		_numParameters = 1;
		_initializeBetas();
	};

	TRecalibrationEMCovariateFunction_intercept(std::vector<std::string> & values){
		_moduleName = RecalModuleFunctionName_intercept;
		_numParameters = 1;
		_initializeBetas();
		_initializValues(values);
	};

	double getEtaTerm(){
		return _betas[0];
	};

	double getEtaTerm(const uint16_t val){
		return _betas[0];
	};
};




//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_polynomial
// A polynomial function
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction_polynomial:public TRecalibrationEMCovariateFunction{
protected:

public:
	TRecalibrationEMCovariateFunction_polynomial(int order){
		_moduleName = RecalModuleFunctionName_polynomial;
		_numParameters = order;
		_initializeBetas();
	};

	TRecalibrationEMCovariateFunction_polynomial(std::vector<std::string> & values){
		_moduleName = RecalModuleFunctionName_polynomial;
		_numParameters = values.size();
		_initializeBetas();
		_initializValues(values);
	};

	double getEtaTerm(const uint16_t val){
		double tmp = _getAsDouble(val);
		return tmp * (_betas[0] + _betas[1] * tmp);
	};
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specific
// A term per discrete value from 0 to maxValue
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction_specific:public TRecalibrationEMCovariateFunction{
protected:
	int _maxValue;

public:
	TRecalibrationEMCovariateFunction_specific(int maxValue){
		_moduleName = RecalModuleFunctionName_specific;
		_maxValue = maxValue;
		_numParameters = _maxValue + 1;
		_initializeBetas();
	};

	TRecalibrationEMCovariateFunction_specific(std::vector<std::string> & values){
		_moduleName = RecalModuleFunctionName_specific;

		//values indicate size
		_numParameters = values.size();
		_maxValue = _numParameters - 1;

		//now initialize betas
		_initializeBetas();

		//now copy values
		_initializValues(values);
	};

	bool checkValueRange(const uint16_t val){
		if(val > _maxValue)
			return false;
		return true;
	};

	double getEtaTerm(const uint16_t val){
		return _betas[val];
	};
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specificMap
// A term per discrete values as indicated with a map
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction_specificMap:public TRecalibrationEMCovariateFunction{
protected:
	int _maxValue;
	uint16_t _mapSize;
	uint16_t* _indexMap;
	bool* _valueUsed;

	void _createIndexMap(){
		_mapSize = _maxValue + 1;
		_indexMap = new uint16_t[_mapSize];
		_valueUsed = new bool[_mapSize];
		for(uint16_t i=0; i<_mapSize; ++i){
			_indexMap[i] = 0;
			_valueUsed[i] = false;
		}
	};

public:
	TRecalibrationEMCovariateFunction_specificMap(std::vector<uint16_t> & values){
		_moduleName = RecalModuleFunctionName_specific;

		//create map based on values
		_numParameters = values.size();
		_initializeBetas();

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

	TRecalibrationEMCovariateFunction_specificMap(std::vector<std::string> & values){
		_moduleName = RecalModuleFunctionName_specific;

		//values indicate size
		_numParameters = values.size();
		_initializeBetas();

		//parse values as pairs separated by a colon (:)
		std::map<uint16_t, double> tmp;
		for(std::string s : values){
			size_t pos = s.find(':');
			if(pos == std::string::npos){
				throw "Can not parse value '" + s + "': missing ':'!";
			}
			tmp.emplace(stringToIntCheck(s.substr(0, pos)), stringToDoubleCheck(s.substr(pos+1)));
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

	bool checkValueRange(const uint16_t val){
		if(!_valueUsed[val])
			return false;
		return true;
	};

	double getEtaTerm(const uint16_t val){
		return _betas[val];
	};
};


#endif /* RECALIBRATION_TRECALIBRATIONEMCOVARIATEFUNCTION_H_ */
