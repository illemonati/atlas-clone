/*
 * TRecalibrationEMModule.h
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#ifndef RECALIBRATION_TRECALIBRATIONEMMODULE_H_
#define RECALIBRATION_TRECALIBRATIONEMMODULE_H_

#include "../stringFunctions.h"
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

//define module names
#define RecalModuleName_none "none"
#define RecalModuleName_intercept "intercept"
#define RecalModuleName_linear "linear"
#define RecalModuleName_quadratic "quadratic"
#define RecalModuleName_specific "specific"


//--------------------------------------------------------------
// TRecalibrationEMTransformationMap
//--------------------------------------------------------------
class TRecalibrationEMTransformationMap{
private:
	uint16_t size;
	uint16_t min, max;
	double* map;
	bool initialized;

public:
	TRecalibrationEMTransformationMap(const uint16_t Min, const uint16_t Max){
		initialized = false;
		min = 0;
		max = 0;
		size = 0;
		map = nullptr;
	};
	TRecalibrationEMTransformationMap(const uint16_t Min, const uint16_t Max){
		min = Min;
		max = Max;
		size = Max - Min + 1;
		map = new double[size];
		for(int i=0; i<size; i++){
			map[i] = 0.0;
		}
	};

	void set(const uint16_t x, const double value){
		map[x - min] = value;
	};

	double operator[](const int x){
		return map[x-min];
	};
};

//--------------------------------------------------------------
// TRecalibrationEMModule
// Base class for recal modules
//--------------------------------------------------------------

class TRecalibrationEMModule{
protected:
	size_t _numParameters;
	std::string _moduleName;
	double* _betas; //betas of the model
	double* _oldBetas; //use during estimation
	bool _initialized;

	//transform values?
	TRecalibrationEMTransformationMap;
	double* transformationMap;

	void _freeBetas();
	void _initializeBetas();
	void _initializValues(std::vector<std::string> & values);

	double _getAsDouble(const uint16_t val){
		if(doTransformation){
			return transformationMap[val];
		} else {
			return (double) val;
		}
	};

	friend class TRecalibrationEMModelNEW;

public:
	TRecalibrationEMModule(){
		_initialized = false;
		_numParameters = 0;
		_betas = nullptr;
		_oldBetas = nullptr;

		doTransformation = false;
		transformationMap = nullptr;
	};

	virtual ~TRecalibrationEMModule(){
		_freeBetas();
	};

	virtual bool checkValueRange(std::vector<std::uint16_t> & values){
		return true;
	};

	void addTransformation(double* pointerToTransformationMap){
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

	virtual bool checkValueRange(uint16_t val){
		return true;
	};
};

//--------------------------------------------------------------
// TRecalibrationEMModule_intercept
// An intercept term
//--------------------------------------------------------------
class TRecalibrationEMModule_intercept:public TRecalibrationEMModule{
protected:

public:
	TRecalibrationEMModule_intercept(){
		_numParameters = 1;
		_initializeBetas();
	};

	TRecalibrationEMModule_intercept(std::vector<std::string> & values){
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
// TRecalibrationEMModule_linear
// A linear function
//--------------------------------------------------------------
class TRecalibrationEMModule_linear:public TRecalibrationEMModule{
protected:

public:
	TRecalibrationEMModule_linear(){
		_numParameters = 1;
		_initializeBetas();
	};

	TRecalibrationEMModule_linear(std::vector<std::string> & values){
		_numParameters = 1;
		_initializeBetas();
		_initializValues(values);
	};

	double getEtaTerm(const uint16_t val){
		return _betas[0] * _getAsDouble(val);
	};
};


//--------------------------------------------------------------
// TRecalibrationEMModule_quadratic
// A quadratic function
//--------------------------------------------------------------
class TRecalibrationEMModule_quadratic:public TRecalibrationEMModule{
protected:

public:
	TRecalibrationEMModule_quadratic(){
		_numParameters = 2;
		_initializeBetas();
	};

	TRecalibrationEMModule_quadratic(std::vector<std::string> & values){
		_numParameters = 2;
		_initializeBetas();
		_initializValues(values);
	};

	double getEtaTerm(const uint16_t val){
		double tmp = _getAsDouble(val);
		return tmp * (_betas[0] + _betas[1] * tmp);
	};
};

//--------------------------------------------------------------
// TRecalibrationEMModule_specific
// A term per discrete value from 0 to maxValue
//--------------------------------------------------------------
class TRecalibrationEMModule_specific:public TRecalibrationEMModule{
protected:
	int _maxValue;

public:
	TRecalibrationEMModule_specific(int maxValue){
		_maxValue = maxValue;
		_numParameters = _maxValue + 1;
		_initializeBetas();
	};

	TRecalibrationEMModule_specific(std::vector<std::string> & values){
		//values indicate size
		_numParameters = values.size();
		_maxValue = _numParameters - 1;

		//now initialize betas
		_initializeBetas();

		//now copy values
		_initializValues(values);
	};

	double getEtaTerm(const uint16_t val){
		return _betas[val];
	};

	virtual bool checkValueRange(uint16_t val){
		if(val < 0 || val > _maxValue)
			return false;
		return true;
	};
};

//--------------------------------------------------------------
// TRecalibrationEMModule_specificMap
// A term per discrete values as indicated with a map
//--------------------------------------------------------------
class TRecalibrationEMModule_specificMap:public TRecalibrationEMModule{
protected:
	int _maxValue;
	MAP!

public:
	TRecalibrationEMModule_specificMap(int maxValue){
		_maxValue = maxValue;
		_numParameters = _maxValue + 1;
		_initializeBetas();
	};

	TRecalibrationEMModule_specificMap(std::vector<std::string> & values){
		//values indicate size
		_numParameters = values.size();
		_maxValue = _numParameters - 1;

		//now initialize betas
		_initializeBetas();

		//now copy values
		_initializValues(values);
	};

	double getEtaTerm(const uint16_t val){
		return _betas[val];
	};

	virtual bool checkValueRange(uint16_t val){
		if(val < 0 || val > _maxValue)
			return false;
		return true;
	};
};


#endif /* RECALIBRATION_TRECALIBRATIONEMMODULE_H_ */
