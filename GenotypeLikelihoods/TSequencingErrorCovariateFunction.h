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
namespace SequencingError {

//--------------------------------------------------------------
// TCovariateFunction
// Base class for recal covariate functions
//--------------------------------------------------------------
class TCovariateFunction{
protected:
	uint16_t _numParameters;
	uint16_t _firstParameterIndex;
	uint16_t _numNonZeroFirstDerivatives;
	uint16_t _numNonZeroSecondDerivatives;
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

	friend class TModel;

public:
	static inline const std::string name = "none";
	TCovariateFunction(){
		_init(0);
	};
	TCovariateFunction(uint16_t FirstParameterIndex){
		_init(FirstParameterIndex);
	};

	virtual ~TCovariateFunction(){};

	void setBeta(uint16_t index, double val){
		if(index < _numParameters)
			_betas[index] = val;
	};

	//size
	uint16_t numParameters() const noexcept { return _numParameters; };
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
	virtual std::string typeString() const noexcept {return name; }
	std::string getModelString() const;
};

//--------------------------------------------------------------
// TCovariateFunction_intercept
// An intercept term
//--------------------------------------------------------------
class TCovariateFunction_intercept:public TCovariateFunction{
protected:
	void _init();

public:
	static inline const std::string name = "intercept";
	TCovariateFunction_intercept();
	TCovariateFunction_intercept(uint16_t FirstParameterIndex);
	TCovariateFunction_intercept(uint16_t FirstParameterIndex, const std::vector<std::string> & values);

	void initialize(uint16_t FirstParameterIndex);
	void initialize(uint16_t FirstParameterIndex, const std::vector<std::string> & values);
	void setIntercept(double val);
	void addToIntercept(double val);
	double getIntercept() const;
	double getEtaTerm() const;
	double getEtaTerm(uint16_t val) const override;
	void fillDerivatives(TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const;
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const override;
	virtual std::string typeString() const noexcept override {return name; }
};

//--------------------------------------------------------------
// TCovariateFunction_polynomial
// A polynomial function
//--------------------------------------------------------------
class TCovariateFunction_polynomial:public TCovariateFunction{
protected:
	void _init(size_t order);

	//TODO: add tmp storage for eta!

public:
	static inline const std::string name = "polynomial";
	TCovariateFunction_polynomial(uint16_t FirstParameterIndex, size_t order);
	TCovariateFunction_polynomial(uint16_t FirstParameterIndex, const std::vector<std::string> & values);

	double getEtaTerm(uint16_t val) const;
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const;
	virtual std::string typeString() const noexcept override {return name; }
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

class TRecalibrationEMCovariateFunction_probit:public TCovariateFunction{
protected:
	uint16_t _secondParameterIndex;
	uint16_t _thirdParameterIndex;

	//tmp storage
	mutable unsigned int _maxValue;
	mutable std::vector<TProbitTmpStorage> _tmpStorage;

	void _init(uint16_t Maxvalue);
	void _expandTmpStorage(uint16_t MaxValue) const;

public:
	static inline const std::string name = "probit";
	TRecalibrationEMCovariateFunction_probit(uint16_t FirstParameterIndex, uint16_t MaxValue);
	TRecalibrationEMCovariateFunction_probit(uint16_t FirstParameterIndex, const std::vector<std::string> & values);

	double getEtaTerm(uint16_t val) const override;
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const override;
	virtual std::string typeString() const noexcept override {return name; }
};

//--------------------------------------------------------------
// TCovariateFunction_specific
// A term per discrete value from 0 to maxValue
//--------------------------------------------------------------
class TCovariateFunction_specific:public TCovariateFunction{
protected:
	uint16_t _maxValue;

	void _init(uint16_t MaxValue);

public:
	static inline const std::string name = "specific";
	TCovariateFunction_specific(uint16_t FirstParameterIndex, uint16_t MaxValue);
	TCovariateFunction_specific(uint16_t FirstParameterIndex, const std::vector<std::string> & betas);

	bool checkValueRange(uint16_t val) const override;
	void adjustValueRanges(const std::vector<uint16_t> & values) override;

	double getEtaTerm(uint16_t val) const override{
		return _betas[val];
	};
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const override;
	double adjustParametersPostEstimation() override { return _normalizeParameters(); };
	virtual std::string typeString() const noexcept override {return name; }
};

//--------------------------------------------------------------
// TCovariateFunction_specificMap
// A term per discrete values as indicated with a map
//--------------------------------------------------------------
struct TCovariateFunctionIndexMapEntry{
	uint16_t index;
	bool used;

	TCovariateFunctionIndexMapEntry(){
		index = 0;
		used = false;
	};
};

class TCovariateFunction_specificMap:public TCovariateFunction{
protected:
	uint16_t _maxValue;
	std::vector<TCovariateFunctionIndexMapEntry> _indexMap; //maps value to parameter index

	void _init(size_t NumParameters);
	void _initMapFromVector(const std::vector<uint16_t> & values);

public:
	static inline const std::string name = "map";
	TCovariateFunction_specificMap(uint16_t FirstParameterIndex, const std::vector<uint16_t> & values);
	TCovariateFunction_specificMap(uint16_t FirstParameterIndex, const std::vector<std::string> & values);
	~TCovariateFunction_specificMap(){};

	bool checkValueRange(uint16_t val) const override;
	void adjustValueRanges(const std::vector<uint16_t> & values) override;

	double getEtaTerm(uint16_t val) const override{
		return _betas[ _indexMap[val].index ];
	};
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second) const override;
	double adjustParametersPostEstimation() override { return _normalizeParameters(); };
	virtual std::string typeString() const noexcept override {return name; }
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_ */
