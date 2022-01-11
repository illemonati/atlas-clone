/*
 * TRecalibrationEMModule.h
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_

#include "stringFunctions.h"
#include "auxiliaryTools.h"
#include "mathFunctions.h"
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

namespace GenotypeLikelihoods{

//define module names
extern const std::string SequencingErrorCovariateFunction_none;
extern const std::string SequencingErrorCovariateFunction_intercept;
extern const std::string SequencingErrorCovariateFunction_polynomial;
extern const std::string SequencingErrorCovariateFunction_probit;
extern const std::string SequencingErrorCovariateFunction_specific;
extern const std::string SequencingErrorCovariateFunction_map;

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

	void _init(uint16_t FirstParameterIndex);
	void _initializeBetas();
	void _initializValues(const std::vector<std::string> & values);
	double _getAsDouble(uint16_t val) const;
	double _normalizeParameters();

	friend class TSequencingErrorModel;

public:
	TSequencingErrorCovariateFunction(){
		_init(0);
	};
	TSequencingErrorCovariateFunction(uint16_t FirstParameterIndex){
		_init(FirstParameterIndex);
	};

	virtual ~TSequencingErrorCovariateFunction(){};

	void setBeta(uint16_t index, const double & val){
		if(index < _numParameters)
			_betas[index] = val;
	};

	//size
	uint16_t numParameters(){ return _numParameters; };
	uint16_t numNonZeroFirstDerivatives(){ return _numNonZeroFirstDerivatives; };
	uint16_t numNonZeroSecondDerivatives(){ return _numNonZeroSecondDerivatives; };

	//check value range: to ensure that data can be recalibrated
	virtual bool checkValueRange(uint16_t ) const { return true; };
	bool checkValueRange(const std::vector<uint16_t> & values) const;

	//adjust value range: to ensure that function can estimated
	virtual void adjustValueRanges(const std::vector<uint16_t> &){};

	//transformation
	void addTransformation(TRecalibrationEMTransformationMap* pointerToTransformationMap){
		doTransformation = true;
		transformationMap = pointerToTransformationMap;
	};

	//estimation
	virtual double getEtaTerm(uint16_t ) const{
		return 0.0;
	};
	void proposeNewParameters(const arma::mat & JxF, uint16_t & index, double & lambda);
	void rejectProposedParameters(){
		for(unsigned int i=0; i<_numParameters; ++i)
			_betas[i] = _oldBetas[i];
	};
	virtual void fillDerivatives(uint16_t , TRecalibrationEMFirstDerivatives &, TRecalibrationEMSecondDerivatives &) const{};
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
	TSequencingErrorCovariateFunction_intercept(uint16_t FirstParameterIndex);
	TSequencingErrorCovariateFunction_intercept(uint16_t FirstParameterIndex, const std::vector<std::string> & values);

	void initialize(uint16_t FirstParameterIndex);
	void initialize(uint16_t FirstParameterIndex, const std::vector<std::string> & values);
	void setIntercept(const double & val);
	void addToIntercept(const double & val);
	double getIntercept() const;
	double getEtaTerm() const;
	double getEtaTerm(uint16_t val) const override;
	void fillDerivatives(TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const;
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const override;
};

//--------------------------------------------------------------
// TSequencingErrorCovariateFunction_polynomial
// A polynomial function
//--------------------------------------------------------------
class TSequencingErrorCovariateFunction_polynomial:public TSequencingErrorCovariateFunction{
protected:
	void _init(const size_t & order);

	//TODO: add tmp storage for eta!

public:
	TSequencingErrorCovariateFunction_polynomial(uint16_t FirstParameterIndex, const size_t & order);
	TSequencingErrorCovariateFunction_polynomial(uint16_t FirstParameterIndex, const std::vector<std::string> & values);

	double getEtaTerm(uint16_t val) const;
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const;
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_probit
// A polynomial function
//--------------------------------------------------------------
class TProbitTmpStorage{
public:
	double _cumulDens_Phi;
	double _normalDens_phi;
	double _eta;
	double _normalDens_q;
	double _normalDens_Beta1;
	double _normalDens_Beta1_q;
	double _normalDens_Beta1_z;
	double _normalDens_Beta1_q_z;
	double _normalDens_Beta1_q2_z;

	TProbitTmpStorage(const std::vector<double> & betas, uint16_t q);
};

class TRecalibrationEMCovariateFunction_probit:public TSequencingErrorCovariateFunction{
protected:
	uint16_t _secondParameterIndex;
	uint16_t _thirdParameterIndex;

	//tmp storage
	mutable unsigned int _maxValue;
	mutable std::vector<TProbitTmpStorage> _tmpStorage;

	void _init(uint16_t Maxvalue);
	void _expandTmpStorage(uint16_t MaxValue) const;

public:
	TRecalibrationEMCovariateFunction_probit(uint16_t FirstParameterIndex, uint16_t MaxValue);
	TRecalibrationEMCovariateFunction_probit(uint16_t FirstParameterIndex, const std::vector<std::string> & values);

	double getEtaTerm(uint16_t val) const override;
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const override;
};

//--------------------------------------------------------------
// TSequencingErrorCovariateFunction_specific
// A term per discrete value from 0 to maxValue
//--------------------------------------------------------------
class TSequencingErrorCovariateFunction_specific:public TSequencingErrorCovariateFunction{
protected:
	uint16_t _maxValue;

	void _init(uint16_t MaxValue);

public:
	TSequencingErrorCovariateFunction_specific(uint16_t FirstParameterIndex, uint16_t MaxValue);
	TSequencingErrorCovariateFunction_specific(uint16_t FirstParameterIndex, const std::vector<std::string> & betas);

	bool checkValueRange(uint16_t val) const override;
	void adjustValueRanges(const std::vector<uint16_t> & values) override;

	double getEtaTerm(uint16_t val) const override{
		return _betas[val];
	};
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const override;
	double adjustParametersPostEstimation() override { return _normalizeParameters(); };
};

//--------------------------------------------------------------
// TSequencingErrorCovariateFunction_specificMap
// A term per discrete values as indicated with a map
//--------------------------------------------------------------
struct TSequencingErrorCovariateFunctionIndexMapEntry{
	uint16_t index;
	bool used;

	TSequencingErrorCovariateFunctionIndexMapEntry(){
		index = 0;
		used = false;
	};
};

class TSequencingErrorCovariateFunction_specificMap:public TSequencingErrorCovariateFunction{
protected:
	uint16_t _maxValue;
	std::vector<TSequencingErrorCovariateFunctionIndexMapEntry> _indexMap; //maps value to parameter index

	void _init(const size_t & NumParameters);
	void _initMapFromVector(const std::vector<uint16_t> & values);

public:
	TSequencingErrorCovariateFunction_specificMap(uint16_t FirstParameterIndex, const std::vector<uint16_t> & values);
	TSequencingErrorCovariateFunction_specificMap(uint16_t FirstParameterIndex, const std::vector<std::string> & values);
	~TSequencingErrorCovariateFunction_specificMap(){};

	bool checkValueRange(uint16_t val) const override;
	void adjustValueRanges(const std::vector<uint16_t> & values) override;

	double getEtaTerm(uint16_t val) const override{
		return _betas[ _indexMap[val].index ];
	};
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const override;
	double adjustParametersPostEstimation() override { return _normalizeParameters(); };
};

}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_ */
