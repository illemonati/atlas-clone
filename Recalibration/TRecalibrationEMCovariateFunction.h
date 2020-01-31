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
	uint16_t _numParameters;
	uint16_t _firstParameterIndex;
	uint16_t _numNonZeroFirstDerivatives;
	uint16_t _numNonZeroSecondDerivatives;
	std::string _moduleName;
	double* _betas; //betas of the model
	double* _oldBetas; //use during estimation
	bool _initialized;

	//transform values?
	TRecalibrationEMTransformationMap* transformationMap;
	bool doTransformation;

	void _init(const uint16_t FirstParameterIndex);
	void _freeBetas();
	void _initializeBetas();
	void _initializValues(std::vector<std::string> & values);
	double _getAsDouble(const uint16_t val);
	double _normalizeParameters();

	friend class TRecalibrationEMModel;

public:
	TRecalibrationEMCovariateFunction(){
		_init(0);
	};
	TRecalibrationEMCovariateFunction(const uint16_t FirstParameterIndex){
		_init(FirstParameterIndex);
	};

	virtual ~TRecalibrationEMCovariateFunction(){ _freeBetas();	};

	void setBeta(uint16_t index, double val){
		if(index < _numParameters)
			_betas[index] = val;
	};
	virtual bool checkValueRange(const uint16_t val){ return true; };
	bool checkValueRange(std::vector<uint16_t> & values);
	uint16_t numParameters(){ return _numParameters; };
	uint16_t numNonZeroFirstDerivatives(){ return _numNonZeroFirstDerivatives; };
	uint16_t numNonZeroSecondDerivatives(){ return _numNonZeroSecondDerivatives; };

	void addTransformation(TRecalibrationEMTransformationMap* pointerToTransformationMap){
		doTransformation = true;
		transformationMap = pointerToTransformationMap;
	};

	virtual double getEtaTerm(const uint16_t val){ return 0.0; }
	void proposeNewParameters(const arma::mat & JxF, uint16_t & index, double & lambda);
	void rejectProposedParameters(){
		for(unsigned int i=0; i<_numParameters; ++i)
			_betas[i] = _oldBetas[i];
	};
	virtual void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second){};
	virtual double adjustParametersPostEstimation(){ return 0.0; };
	std::string getModelString();
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_intercept
// An intercept term
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction_intercept:public TRecalibrationEMCovariateFunction{
protected:
	void _init();

public:
	TRecalibrationEMCovariateFunction_intercept();
	TRecalibrationEMCovariateFunction_intercept(const uint16_t FirstParameterIndex);
	TRecalibrationEMCovariateFunction_intercept(const uint16_t FirstParameterIndex, std::vector<std::string> & values);

	void addToIntercept(const double val){
		_betas[0] += val;
	};

	double getEtaTerm(){
		return _betas[0];
	};

	double getEtaTerm(const uint16_t val){
		return _betas[0];
	};

	void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second);
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_polynomial
// A polynomial function
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction_polynomial:public TRecalibrationEMCovariateFunction{
protected:
	void _init(const size_t order);

public:
	TRecalibrationEMCovariateFunction_polynomial(const uint16_t FirstParameterIndex, const size_t order);
	TRecalibrationEMCovariateFunction_polynomial(const uint16_t FirstParameterIndex, std::vector<std::string> & values);

	double getEtaTerm(const uint16_t val);
	void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second);
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specific
// A term per discrete value from 0 to maxValue
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction_specific:public TRecalibrationEMCovariateFunction{
protected:
	unsigned int _maxValue;

	void _init(const size_t MaxValue);

public:
	TRecalibrationEMCovariateFunction_specific(const uint16_t FirstParameterIndex, size_t MaxValue);
	TRecalibrationEMCovariateFunction_specific(const uint16_t FirstParameterIndex, std::vector<std::string> & values);

	bool checkValueRange(const uint16_t val){
		if(val > _maxValue)
			return false;
		return true;
	};

	double getEtaTerm(const uint16_t val){
		return _betas[val];
	};
	void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second);
	double adjustParametersPostEstimation(){ return _normalizeParameters(); };
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_specificMap
// A term per discrete values as indicated with a map
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction_specificMap:public TRecalibrationEMCovariateFunction{
protected:
	unsigned int _maxValue;
	uint16_t _mapSize;
	uint16_t* _indexMap; //maps value to parameter index
	bool* _valueUsed;

	void _init(size_t NumParameters);
	void _createIndexMap();
	void _freeIndexMap();

public:
	TRecalibrationEMCovariateFunction_specificMap(const uint16_t FirstParameterIndex, std::vector<uint16_t> & values);
	TRecalibrationEMCovariateFunction_specificMap(const uint16_t FirstParameterIndex, std::vector<std::string> & values);
	~TRecalibrationEMCovariateFunction_specificMap(){ _freeIndexMap(); };

	bool checkValueRange(const uint16_t val){
		if(!_valueUsed[val])
			return false;
		return true;
	};

	double getEtaTerm(const uint16_t val){
		return _betas[ _indexMap[val] ];
	};
	void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second);
	double adjustParametersPostEstimation(){ return _normalizeParameters(); };
};


#endif /* RECALIBRATION_TRECALIBRATIONEMCOVARIATEFUNCTION_H_ */
