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


//--------------------------------------------------------------------
// Derivatives
//--------------------------------------------------------------------
struct TRecalibrationEMFirstDerivative{
	size_t index;
	double derivative;
};

typedef std::vector<TRecalibrationEMFirstDerivative> TRecalibrationEMFirstDerivatives;

struct TRecalibrationEMSecondDerivative{
	size_t index1;
	size_t index2;
	double derivative;
};

typedef std::vector<TRecalibrationEMSecondDerivative> TRecalibrationEMSecondDerivatives;

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction
// Base class for recal covariate functions
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction{
protected:
	size_t _numParameters;
	size_t _firstParameterIndex;
	size_t _numNonZeroFirstDerivatives;
	size_t _numNonZeroSecondDerivatives;
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
	double _getAsDouble(const uint16_t val);

	friend class TRecalibrationEMModel;

public:
	TRecalibrationEMCovariateFunction(const size_t FirstParameterIndex){
		_moduleName = RecalModuleFunctionName_none;
		_initialized = false;
		_numParameters = 0;
		_firstParameterIndex = FirstParameterIndex;
		_numNonZeroFirstDerivatives = 0;
		_numNonZeroSecondDerivatives = 0;
		_betas = nullptr;
		_oldBetas = nullptr;

		doTransformation = false;
		transformationMap = nullptr;
	};

	virtual ~TRecalibrationEMCovariateFunction(){ _freeBetas();	};
	virtual bool checkValueRange(const uint16_t val){ return true; };
	bool checkValueRange(std::vector<uint16_t> & values);
	int numParameters(){ return _numParameters; };
	int numNonZeroFirstDerivatives(){ return _numNonZeroFirstDerivatives; };
	int numNonZeroSecondDerivatives(){ return _numNonZeroSecondDerivatives; };

	void addTransformation(TRecalibrationEMTransformationMap* pointerToTransformationMap){
		doTransformation = true;
		transformationMap = pointerToTransformationMap;
	};

	virtual double getEtaTerm(const uint16_t val){ return 0.0; }
	void proposeNewParameters(const arma::mat & JxF, size_t & index, double & lambda);
	void rejectProposedParameters(){
		for(unsigned int i=0; i<_numParameters; ++i)
			_betas[i] = _oldBetas[i];
	};
	virtual void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, size_t & indexFirst, TRecalibrationEMSecondDerivatives & second, size_t & indexsecond){};
	std::string getModelString();
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_intercept
// An intercept term
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction_intercept:public TRecalibrationEMCovariateFunction{
protected:
	void _init(){
		_moduleName = RecalModuleFunctionName_intercept;
		_numParameters = 1;
		_numNonZeroFirstDerivatives = 1;
		_numNonZeroSecondDerivatives = 0;
		_initializeBetas();
	};

public:
	TRecalibrationEMCovariateFunction_intercept(const size_t FirstParameterIndex){
		_init();
	};

	TRecalibrationEMCovariateFunction_intercept(const size_t FirstParameterIndex, std::vector<std::string> & values){
		_moduleName = RecalModuleFunctionName_intercept;
		_init();
		_initializValues(values);
	};

	double getEtaTerm(){
		return _betas[0];
	};

	double getEtaTerm(const uint16_t val){
		return _betas[0];
	};

	void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, size_t & indexFirst, TRecalibrationEMSecondDerivatives & second, size_t & indexsecond);
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_polynomial
// A polynomial function
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction_polynomial:public TRecalibrationEMCovariateFunction{
protected:
	void _init(const size_t order){
		_moduleName = RecalModuleFunctionName_polynomial;
		_numParameters = order;
		_numNonZeroFirstDerivatives = order;
		_numNonZeroSecondDerivatives = 0;
		_initializeBetas();
	};

public:
	TRecalibrationEMCovariateFunction_polynomial(const size_t FirstParameterIndex, int order){
		_init(order);
	};

	TRecalibrationEMCovariateFunction_polynomial(const size_t FirstParameterIndex, std::vector<std::string> & values){
		_init(values.size());
		_initializValues(values);
	};

	double getEtaTerm(const uint16_t val){
		double tmp = _getAsDouble(val);
		return tmp * (_betas[0] + _betas[1] * tmp);
	};

	void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, size_t & indexFirst, TRecalibrationEMSecondDerivatives & second, size_t & indexsecond);
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specific
// A term per discrete value from 0 to maxValue
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction_specific:public TRecalibrationEMCovariateFunction{
protected:
	size_t _maxValue;

	void _init(const size_t MaxValue);

public:
	TRecalibrationEMCovariateFunction_specific(const size_t FirstParameterIndex, size_t MaxValue);
	TRecalibrationEMCovariateFunction_specific(const size_t FirstParameterIndex, std::vector<std::string> & values);

	bool checkValueRange(const uint16_t val){
		if(val > _maxValue)
			return false;
		return true;
	};

	double getEtaTerm(const uint16_t val){
		return _betas[val];
	};
	void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, size_t & indexFirst, TRecalibrationEMSecondDerivatives & second, size_t & indexsecond);
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specificMap
// A term per discrete values as indicated with a map
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction_specificMap:public TRecalibrationEMCovariateFunction{
protected:
	int _maxValue;
	uint16_t _mapSize;
	uint16_t* _indexMap; //maps value to parameter index
	bool* _valueUsed;

	void _init(size_t NumParameters);
	void _createIndexMap();
	void _freeIndexMap();

public:
	TRecalibrationEMCovariateFunction_specificMap(const size_t FirstParameterIndex, std::vector<uint16_t> & values);
	TRecalibrationEMCovariateFunction_specificMap(const size_t FirstParameterIndex, std::vector<std::string> & values);
	~TRecalibrationEMCovariateFunction_specificMap(){ _freeIndexMap(); };
	bool checkValueRange(const uint16_t val){
		if(!_valueUsed[val])
			return false;
		return true;
	};

	double getEtaTerm(const uint16_t val){
		return _betas[val];
	};
	void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, size_t & indexFirst, TRecalibrationEMSecondDerivatives & second, size_t & indexsecond);
};


#endif /* RECALIBRATION_TRECALIBRATIONEMCOVARIATEFUNCTION_H_ */
