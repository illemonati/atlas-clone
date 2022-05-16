/*
 * TRecalibrationEMModule.h
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_

#include <stddef.h>
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

#include "auxiliaryTools.h"

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
	TFunction(uint16_t FirstParameterIndex = 0) : _firstParameterIndex(FirstParameterIndex) {}
	virtual ~TFunction() = default;

	// non-virtuals
	void proposeNewParameters(const arma::mat &JxF, uint16_t &index, double &lambda) noexcept;
	void rejectProposedParameters() noexcept;
	constexpr uint16_t firstParameterIndex() const noexcept { return _firstParameterIndex; };

	// virtuals
	virtual double &beta(uint16_t index) noexcept      = 0;
	virtual double beta(uint16_t index) const noexcept = 0;

	virtual bool recalibrates() const noexcept = 0;

	virtual uint16_t numParameters() const noexcept               = 0;
	virtual uint16_t numNonZeroFirstDerivatives() const noexcept  = 0;
	virtual uint16_t numNonZeroSecondDerivatives() const noexcept = 0;

	// check value range: to ensure that data can be recalibrated
	virtual bool checkOrInitValueRange(const std::vector<uint16_t> &values) = 0;

	// estimation
	virtual double getEtaTerm(uint16_t val) const noexcept                                 = 0;
	virtual void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives &first,
								 TRecalibrationEMSecondDerivatives &second) const noexcept = 0;
	virtual double adjustParametersPostEstimation() noexcept                               = 0;
	virtual std::string typeString() const noexcept                                        = 0;
	std::string getModelString() const;
};

//--------------------------------------------------------------
// TCovariateFunction_intercept
// An intercept term
//--------------------------------------------------------------
class TIntercept : public TFunction {
private:
	double _beta    = 0.;
	double _oldbeta = 0.;

protected:
	double &_oldBeta(uint16_t) noexcept override { return _oldbeta; }
	double _oldBeta(uint16_t) const noexcept override { return _oldbeta; }
public:
	static inline const std::string name = "intercept";
	// No beta
	TIntercept(uint16_t FirstParameterIndex): TFunction(FirstParameterIndex) {};
	// With beta
	TIntercept(uint16_t FirstParameterIndex, const std::string &beta);

	bool recalibrates() const noexcept override {return _beta != 0.;}

	uint16_t numParameters() const noexcept override { return 1; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 1; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	double &beta(uint16_t) noexcept override { return _beta; }
	double beta(uint16_t) const noexcept override { return _beta; }

	bool checkOrInitValueRange(const std::vector<uint16_t>&) noexcept override {return true;};

	constexpr double intercept() const noexcept { return _beta; }
	constexpr double &intercept() noexcept { return _beta; }

	double getEtaTerm(uint16_t = 0) const noexcept override { return _beta; };
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives &first,
			     TRecalibrationEMSecondDerivatives &second) const noexcept override;
	double adjustParametersPostEstimation() noexcept override { return 0.; }

	virtual std::string typeString() const noexcept override;
};

//--------------------------------------------------------------
// TCovariateFunction_polynomial
// A polynomial function
//--------------------------------------------------------------
template <size_t O>
class TPolynomial : public TFunction {
	static_assert(O > 0);
private:
	std::array<double, O> _betas;    // betas of the model
	std::array<double, O> _oldBetas; // use during estimation
	bool _recal;
	TRecalibrationEMTransformationMap *_transformationMap = nullptr;

	double _getAsDouble(uint16_t val) const noexcept {
		return _transformationMap ? _transformationMap->get(val) : val;
	}

protected:
	double &_oldBeta(uint16_t i) noexcept override { return _oldBetas[i]; }
	double _oldBeta(uint16_t i) const noexcept override { return _oldBetas[i]; }
public:
	static inline const std::string name = "polynomial";

	// No betas
	TPolynomial(uint16_t FirstParameterIndex, TRecalibrationEMTransformationMap *transformationMap = nullptr)
		: TFunction(FirstParameterIndex), _recal(false), _transformationMap(transformationMap) {
		_betas.front() = 1.;
	}
	// With betas
	TPolynomial(uint16_t FirstParameterIndex, const std::vector<std::string> &betas,
				TRecalibrationEMTransformationMap *transformationMap = nullptr)
		: TFunction(FirstParameterIndex), _recal(true), _transformationMap(transformationMap) {
		std::transform(betas.cbegin(), betas.cend(), _betas.begin(), coretools::str::convertStringCheck<double>);
	}

	bool recalibrates() const noexcept override {return _recal;}

	uint16_t numParameters() const noexcept override { return O; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return O; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	double &beta(uint16_t i) noexcept override { return _betas[i]; }
	double beta(uint16_t i) const noexcept override { return _betas[i]; }

	bool checkOrInitValueRange(const std::vector<uint16_t> &values) noexcept override {
		return !_transformationMap ||
			   std::all_of(values.begin(), values.end(), [tm = this->_transformationMap](auto v) { return tm->checkRange(v); });
	};

	double adjustParametersPostEstimation() noexcept override { return 0.; }

	double getEtaTerm(uint16_t val) const noexcept override {
		const double v = _getAsDouble(val);
		if constexpr (O == 1) {
			return _betas.front() * v;
		} else if constexpr (O == 2) {
			return v * (_betas[0] + v * _betas[1]);
		} else if constexpr (O == 3) {
			return v * (_betas[0] + v * (_betas[1] + _betas[2] * v));
		} else {
			double vpi = v;
			double sum = _betas[0] * vpi;

			for (size_t i = 1; i < O; ++i) {
				vpi *= v;
				sum += _betas[i] * vpi;
			}
			return sum;
		}
	};

	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives &first,
						 TRecalibrationEMSecondDerivatives &) const noexcept override {
		const double v = _getAsDouble(val);

		if constexpr (O == 1) {
			first.add(firstParameterIndex(), v);
		} else if constexpr (O == 2) {
			first.add(firstParameterIndex(), v);
			first.add(firstParameterIndex() + 1, v*v);
		} else if constexpr (O == 3) {
			first.add(firstParameterIndex(), v);
			first.add(firstParameterIndex() + 1, v*v);
			first.add(firstParameterIndex() + 2, v*v*v);
		} else {
			double vpi = v;
			first.add(firstParameterIndex(), vpi);
			for (size_t i = 1; i < numParameters(); ++i) {
				vpi *= v;
				first.add(firstParameterIndex() + i, vpi);
			}
		}
	}
	std::string typeString() const noexcept override {
		using coretools::str::toString;
		std::string s = name + '(' + toString(O) + ')';
		s += '[';
		for (const auto b : _betas) { s += toString(b) + ','; }
		return s.substr(0, s.size() - 1) + ']';
	}
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
	bool _recal;

	// tmp storage
	mutable std::vector<TProbitTmpStorage> _tmpStorage;

	void _expandTmpStorage(uint16_t MaxValue) const;

protected:
	double &_oldBeta(uint16_t i) noexcept override { return _oldBetas[i]; }
	double _oldBeta(uint16_t i) const noexcept override { return _oldBetas[i]; }

public:
	static inline const std::string name = "probit";
	// No betas
	TProbit(uint16_t FirstParameterIndex);
	// With betas
	TProbit(uint16_t FirstParameterIndex, const std::vector<std::string> &betas);

	bool recalibrates() const noexcept override {return _recal;}

	uint16_t numParameters() const noexcept override { return 3; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 3; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 6; };

	bool checkOrInitValueRange(const std::vector<uint16_t>&) noexcept override {return true;};

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
	std::vector<double> _betas;    // betas of the model
	std::vector<double> _oldBetas; // use during estimation
	bool _recal;

	void _init(uint16_t MaxValue);
	bool _checkValueRange(uint16_t val) const noexcept { return val < numParameters(); };
	void _adjustValueRanges(const std::vector<uint16_t> &values);
protected:
	double &_oldBeta(uint16_t i) noexcept override { return _oldBetas[i]; }
	double _oldBeta(uint16_t i) const noexcept override { return _oldBetas[i]; }
public:
	static inline const std::string name = "specific";

	// No Range nor betas
	TSpecific(uint16_t FirstParameterIndex);
	// Range and betas
	TSpecific(uint16_t FirstParameterIndex, const std::vector<std::string> &betas);

	bool recalibrates() const noexcept override {return _recal;}

	uint16_t numParameters() const noexcept override { return _betas.size(); };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 1; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	double &beta(uint16_t i) noexcept override { return _betas[i]; }
	double beta(uint16_t i) const noexcept override { return _betas[i]; }

	bool checkOrInitValueRange(const std::vector<uint16_t> &values) override {
		if (numParameters() == 0) {
			_adjustValueRanges(values);
			return true;
		}

		for (uint16_t val : values) {
			if (!_checkValueRange(val)) return false;
		}
		return true;
	}

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
	std::vector<TIndexMapEntry> _indexMap; // maps value to parameter index
	bool _recal;

	void _init(size_t NumParameters);
	void _initMapFromVector(const std::vector<uint16_t> &values);
	bool _checkValueRange(uint16_t val) const noexcept { return val < _indexMap.size() && _indexMap[val].used; };
	void _adjustValueRanges(const std::vector<uint16_t> &values);
protected:
	double &_oldBeta(uint16_t i) noexcept override { return _oldBetas[i]; }
	double _oldBeta(uint16_t i) const noexcept override { return _oldBetas[i]; }

public:
	static inline const std::string name = "map";
	// No Range nor betas
	TSpecificMap(uint16_t FirstParameterIndex);
	// betas and range
	TSpecificMap(uint16_t FirstParameterIndex, const std::vector<std::string> &betas);

	bool recalibrates() const noexcept override {return _recal;}

	uint16_t numParameters() const noexcept override { return _betas.size(); };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 1; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	double &beta(uint16_t i) noexcept override { return _betas[i]; }
	double beta(uint16_t i) const noexcept override { return _betas[i]; }

	bool checkOrInitValueRange(const std::vector<uint16_t> &values) override {
		if (numParameters() == 0) {
			_adjustValueRanges(values);
			return true;
		}
		for (uint16_t val : values) {
			if (!_checkValueRange(val)) return false;
		}
		return true;
	}

	double adjustParametersPostEstimation() noexcept override;

	double getEtaTerm(uint16_t val) const noexcept override { return _betas[_indexMap[val].index]; };
	void fillDerivatives(uint16_t val, TRecalibrationEMFirstDerivatives &first,
			     TRecalibrationEMSecondDerivatives &second) const noexcept override;
	virtual std::string typeString() const noexcept override { return name; }
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_ */
