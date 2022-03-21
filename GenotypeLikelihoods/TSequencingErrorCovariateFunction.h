/*
 * TRecalibrationEMModule.h
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_

#include "auxiliaryTools.h"
#include "mathFunctions.h"
#include "stringFunctions.h"
#include <cstdint>
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------
// TCovariateFunction
// Base class for recal covariate functions
//--------------------------------------------------------------
class TFunction {
private:
	uint16_t _firstParameterIndex;
protected:
	virtual double &_oldBeta(uint16_t index) noexcept      = 0;
	virtual double _oldBeta(uint16_t index) const noexcept = 0;
public:
	static inline const std::string name = "none";

	TFunction(uint16_t FirstParameterIndex = 0) : _firstParameterIndex(FirstParameterIndex) {}
	virtual ~TFunction()                               = default;

	virtual double &beta(uint16_t index) noexcept      = 0;
	virtual double beta(uint16_t index) const noexcept = 0;

	// size
	constexpr uint16_t firstParameterIndex() const noexcept { return _firstParameterIndex; };
	virtual uint16_t numParameters() const noexcept               = 0;
	virtual uint16_t numNonZeroFirstDerivatives() const noexcept  = 0;
	virtual uint16_t numNonZeroSecondDerivatives() const noexcept = 0;

	// check value range: to ensure that data can be recalibrated
	virtual bool checkValueRange(uint16_t) const noexcept = 0;
	bool checkValueRange(const std::vector<uint16_t> &values) const noexcept;

	// estimation
	virtual double getEtaTerm(uint16_t) const noexcept = 0;
	void proposeNewParameters(const arma::mat &JxF, uint16_t &index, double &lambda) noexcept;
	void rejectProposedParameters() noexcept;
	virtual void fillDerivatives(uint16_t, TRecalibrationEMFirstDerivatives &,
				     TRecalibrationEMSecondDerivatives &) const noexcept = 0;
	virtual double adjustParametersPostEstimation() noexcept                         = 0;
	virtual std::string typeString() const noexcept                                  = 0;
	std::string getModelString() const;
};

//--------------------------------------------------------------
// TCovariateFunction_intercept
// An intercept term
//--------------------------------------------------------------
class TIntercept : public TFunction {
private:
	double _beta;
	double _oldbeta;

protected:
	double &_oldBeta(uint16_t) noexcept override { return _oldbeta; }
	double _oldBeta(uint16_t) const noexcept override { return _oldbeta; }

public:
	static inline const std::string name = "intercept";
	TIntercept(uint16_t FirstParameterIndex, const std::string &value);

	uint16_t numParameters() const noexcept override { return 1; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 1; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	double &beta(uint16_t) noexcept override { return _beta; }
	double beta(uint16_t) const noexcept override { return _beta; }

	bool checkValueRange(uint16_t) const noexcept override { return true; }
	constexpr double intercept() const noexcept { return _beta; }
	constexpr double &intercept() noexcept { return _beta; }


	double getEtaTerm(uint16_t = 0) const noexcept override { return _beta; };
	void fillDerivatives(TRecalibrationEMFirstDerivatives &first,
			     TRecalibrationEMSecondDerivatives &second) const noexcept;
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives &first,
			     TRecalibrationEMSecondDerivatives &second) const noexcept override;
	double adjustParametersPostEstimation() noexcept override { return 0.; }

	virtual std::string typeString() const noexcept override { return name; }
};

//--------------------------------------------------------------
// TCovariateFunction_polynomial
// A polynomial function
//--------------------------------------------------------------
class TPolynomial : public TFunction {
private:
	std::vector<double> _betas;    // betas of the model
	std::vector<double> _oldBetas; // use during estimation
	void _init(size_t order);
	uint16_t _order;
	TRecalibrationEMTransformationMap *_transformationMap = nullptr;
	double _getAsDouble(uint16_t val) const noexcept {
		return _transformationMap ? _transformationMap->get(val) : val;
	};

protected:
	double &_oldBeta(uint16_t i) noexcept override { return _oldBetas[i]; }
	double _oldBeta(uint16_t i) const noexcept override { return _oldBetas[i]; }
	// TODO: add tmp storage for eta!

public:
	static inline const std::string name = "polynomial";
	TPolynomial(uint16_t FirstParameterIndex, size_t order,
		    TRecalibrationEMTransformationMap *transformationMap = nullptr);
	TPolynomial(uint16_t FirstParameterIndex, const std::vector<std::string> &values,
		    TRecalibrationEMTransformationMap *transformationMap = nullptr);

	uint16_t numParameters() const noexcept override { return _order; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return _order; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	double &beta(uint16_t i) noexcept override { return _betas[i]; }
	double beta(uint16_t i) const noexcept override { return _betas[i]; }


	bool checkValueRange(uint16_t) const noexcept override;
	double adjustParametersPostEstimation() noexcept override { return 0.; }

	double getEtaTerm(uint16_t val) const noexcept override;
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives &first,
			     TRecalibrationEMSecondDerivatives &second) const noexcept override;
	virtual std::string typeString() const noexcept override { return name; }
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_probit
// A polynomial function
//--------------------------------------------------------------

class TProbit : public TFunction {
private:
	struct TProbitTmpStorage {
		double cumulDens_Phi;
		double normalDens_phi;
		double eta;
		double normalDens_q;
		double normalDens_Beta1;
		double normalDens_Beta1_q;
		double normalDens_Beta1_z;
		double normalDens_Beta1_q_z;
		double normalDens_Beta1_q2_z;

		TProbitTmpStorage(const std::array<double, 3> &betas, uint16_t q);
	};
	std::array<double, 3> _betas;    // betas of the model
	std::array<double, 3> _oldBetas; // use during estimation

	// tmp storage
	mutable unsigned int _maxValue;
	mutable std::vector<TProbitTmpStorage> _tmpStorage;

	void _init(uint16_t Maxvalue);
	void _expandTmpStorage(uint16_t MaxValue) const;

protected:
	double &_oldBeta(uint16_t i) noexcept override { return _oldBetas[i]; }
	double _oldBeta(uint16_t i) const noexcept override { return _oldBetas[i]; }

public:
	static inline const std::string name = "probit";
	TProbit(uint16_t FirstParameterIndex, uint16_t MaxValue);
	TProbit(uint16_t FirstParameterIndex, const std::vector<std::string> &values);

	uint16_t numParameters() const noexcept override { return 3; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 3; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 6; };

	double &beta(uint16_t i) noexcept override { return _betas[i]; }
	double beta(uint16_t i) const noexcept override { return _betas[i]; }

	double getEtaTerm(uint16_t val) const noexcept override;
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives &first,
			     TRecalibrationEMSecondDerivatives &second) const noexcept override;
	virtual std::string typeString() const noexcept override { return name; }
};

//--------------------------------------------------------------
// TCovariateFunction_specific
// A term per discrete value from 0 to maxValue
//--------------------------------------------------------------
class TSpecific : public TFunction {
private:
	uint16_t _maxValue;
	std::vector<double> _betas;    // betas of the model
	std::vector<double> _oldBetas; // use during estimation

	void _init(uint16_t MaxValue);

protected:
	double &_oldBeta(uint16_t i) noexcept override { return _oldBetas[i]; }
	double _oldBeta(uint16_t i) const noexcept override { return _oldBetas[i]; }

public:
	static inline const std::string name = "specific";
	TSpecific(uint16_t FirstParameterIndex, uint16_t MaxValue);
	TSpecific(uint16_t FirstParameterIndex, const std::vector<std::string> &betas);

	uint16_t numParameters() const noexcept override { return _maxValue + 1; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 1; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	double &beta(uint16_t i) noexcept override { return _betas[i]; }
	double beta(uint16_t i) const noexcept override { return _betas[i]; }

	bool checkValueRange(uint16_t val) const noexcept override { return val <= _maxValue; };
	void adjustValueRanges(const std::vector<uint16_t> &values);

	double adjustParametersPostEstimation() noexcept override;

	double getEtaTerm(uint16_t val) const noexcept override { return _betas[val]; };
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives &first,
			     TRecalibrationEMSecondDerivatives &second) const noexcept override;
	virtual std::string typeString() const noexcept override { return name; }
};

//--------------------------------------------------------------
// TCovariateFunction_specificMap
// A term per discrete values as indicated with a map
//--------------------------------------------------------------
struct TIndexMapEntry {
	uint16_t index = 0;
	bool used      = false;
};

class TSpecificMap : public TFunction {
private:
	std::vector<double> _betas;    // betas of the model
	std::vector<double> _oldBetas; // use during estimation
	uint16_t _numParameters;
	uint16_t _maxValue;
	std::vector<TIndexMapEntry> _indexMap; // maps value to parameter index

	void _init(size_t NumParameters);
	void _initMapFromVector(const std::vector<uint16_t> &values);

protected:
	double &_oldBeta(uint16_t i) noexcept override { return _oldBetas[i]; }
	double _oldBeta(uint16_t i) const noexcept override { return _oldBetas[i]; }

public:
	static inline const std::string name = "map";
	TSpecificMap(uint16_t FirstParameterIndex, const std::vector<uint16_t> &values);
	TSpecificMap(uint16_t FirstParameterIndex, const std::vector<std::string> &values);
	~TSpecificMap(){};

	uint16_t numParameters() const noexcept override { return _numParameters; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 1; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	double &beta(uint16_t i) noexcept override { return _betas[i]; }
	double beta(uint16_t i) const noexcept override { return _betas[i]; }

	bool checkValueRange(uint16_t val) const noexcept override { return val <= _maxValue || _indexMap[val].used; };
	void adjustValueRanges(const std::vector<uint16_t> &values);

	double adjustParametersPostEstimation() noexcept override;

	double getEtaTerm(uint16_t val) const noexcept override { return _betas[_indexMap[val].index]; };
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives &first,
			     TRecalibrationEMSecondDerivatives &second) const noexcept override;
	virtual std::string typeString() const noexcept override { return name; }
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_ */
