/*
 * TRecalibrationEMModule.h
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_

#include "stringFunctions.h"

#include "../TNormalDistribution.h"
#include "auxiliaryTools.h"
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>


//define module names
#define SequencingErrorCovariateFunction_none "none"
#define SequencingErrorCovariateFunction_intercept "intercept"
#define SequencingErrorCovariateFunction_polynomial "polynomial"
#define SequencingErrorCovariateFunction_probit "probit"
#define SequencingErrorCovariateFunction_specific "specific"

namespace GenotypeLikelihoods{

//--------------------------------------------------------------
// TSequencingErrorCovariateFunction
// Base class for recal covariate functions
//--------------------------------------------------------------
class TSequencingErrorCovariateFunction{
protected:
	uint16_t _numParameters;
	uint16_t _firstParameterIndex;
	uint16_t _numNonZeroFirstDerivatives;
	uint16_t _numNonZeroSecondDerivatives;
	std::string _moduleName;
	std::vector<double> _betas; //betas of the model
	std::vector<double> _oldBetas; //use during estimation

	//transform values?
	TRecalibrationEMTransformationMap* transformationMap;
	bool doTransformation;

	void _init(const uint16_t FirstParameterIndex);
	void _initializeBetas();
	void _initializValues(std::vector<std::string> & values);
	double _getAsDouble(const uint16_t val) const;
	double _normalizeParameters();

	friend class TSequencingErrorModel;

public:
	TSequencingErrorCovariateFunction(){
		_init(0);
	};
	TSequencingErrorCovariateFunction(const uint16_t FirstParameterIndex){
		_init(FirstParameterIndex);
	};

	virtual ~TSequencingErrorCovariateFunction(){};

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

	virtual double getEtaTerm(const uint16_t val){
		return 0.0;
	};
	void proposeNewParameters(const arma::mat & JxF, uint16_t & index, double & lambda);
	void rejectProposedParameters(){
		for(unsigned int i=0; i<_numParameters; ++i)
			_betas[i] = _oldBetas[i];
	};
	virtual void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second){};
	virtual double adjustParametersPostEstimation(){ return 0.0; };
	std::string getModelString() const;
};

//--------------------------------------------------------------
// TSequencingErrorCovariateFunction_intercept
// An intercept term
//--------------------------------------------------------------
class TSequencingErrorCovariateFunction_intercept:public TSequencingErrorCovariateFunction{
protected:
	void _init();

public:
	TSequencingErrorCovariateFunction_intercept();
	TSequencingErrorCovariateFunction_intercept(const uint16_t FirstParameterIndex);
	TSequencingErrorCovariateFunction_intercept(const uint16_t FirstParameterIndex, std::vector<std::string> & values);

	void initialize(const uint16_t FirstParameterIndex, std::vector<std::string> & values);
	void setIntercept(const double val);
	void addToIntercept(const double val);
	double getIntercept() const;
	double getEtaTerm();
	double getEtaTerm(const uint16_t val);
	void fillDerivatives(TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second);
	void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second);
};

//--------------------------------------------------------------
// TSequencingErrorCovariateFunction_polynomial
// A polynomial function
//--------------------------------------------------------------
class TSequencingErrorCovariateFunction_polynomial:public TSequencingErrorCovariateFunction{
protected:
	void _init(const size_t order);

	//TODO: add tmp storage for eta!

public:
	TSequencingErrorCovariateFunction_polynomial(const uint16_t FirstParameterIndex, const size_t order);
	TSequencingErrorCovariateFunction_polynomial(const uint16_t FirstParameterIndex, std::vector<std::string> & values);

	double getEtaTerm(const uint16_t val);
	void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second);
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_probit
// A polynomial function
//--------------------------------------------------------------
class TRecalibrationEMCovariateFunction_probit:public TSequencingErrorCovariateFunction{
protected:
	unsigned int _maxValue;
	TNormalDistribution _normalDist;
	uint16_t _secondParameterIndex;
	uint16_t _thirdParameterIndex;

	//tmp storage
	//TODO: make for all thing needed!
	double* _cumulDens_Phi;
	double* _normalDens_phi;
	double* _eta;
	double* _normalDens_q;
	double* _normalDens_Beta1;
	double* _normalDens_Beta1_q;
	double* _normalDens_Beta1_z;
	double* _normalDens_Beta1_q_z;
	double* _normalDens_Beta1_q2_z;

	bool _tmpStorageInitialized;

	void _init(const uint16_t Maxvalue);
	void _allocateTmpStorage();
	void _freeTmpStorage();
	void _fillTmpStorage();
	void _expandTmpStorage(const uint16_t MaxValue);

public:
	TRecalibrationEMCovariateFunction_probit(const uint16_t FirstParameterIndex, uint16_t MaxValue);
	TRecalibrationEMCovariateFunction_probit(const uint16_t FirstParameterIndex, std::vector<std::string> & values);
	~TRecalibrationEMCovariateFunction_probit(){ _freeTmpStorage(); };

	double getEtaTerm(const uint16_t val);
	void fillDerivatives(const uint16_t & val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second);
};


//--------------------------------------------------------------
// TSequencingErrorCovariateFunction_specific
// A term per discrete value from 0 to maxValue
//--------------------------------------------------------------
class TSequencingErrorCovariateFunction_specific:public TSequencingErrorCovariateFunction{
protected:
	uint16_t _maxValue;

	void _init(const uint16_t MaxValue);

public:
	TSequencingErrorCovariateFunction_specific(const uint16_t FirstParameterIndex, const uint16_t MaxValue);
	TSequencingErrorCovariateFunction_specific(const uint16_t FirstParameterIndex, std::vector<std::string> & values);

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
// TSequencingErrorCovariateFunction_specificMap
// A term per discrete values as indicated with a map
//--------------------------------------------------------------
class TSequencingErrorCovariateFunction_specificMap:public TSequencingErrorCovariateFunction{
protected:
	uint16_t _maxValue;
	uint16_t _mapSize;
	uint16_t* _indexMap; //maps value to parameter index
	bool* _valueUsed;

	void _init(size_t NumParameters);
	void _createIndexMap();
	void _freeIndexMap();

public:
	TSequencingErrorCovariateFunction_specificMap(const uint16_t FirstParameterIndex, std::vector<uint16_t> & values);
	TSequencingErrorCovariateFunction_specificMap(const uint16_t FirstParameterIndex, std::vector<std::string> & values);
	~TSequencingErrorCovariateFunction_specificMap(){ _freeIndexMap(); };

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

}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_ */
