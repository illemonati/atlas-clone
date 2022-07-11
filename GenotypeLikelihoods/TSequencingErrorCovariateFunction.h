/*
 * TRecalibrationEMModule.h
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_

#include <array>
#include <cstdint>
#include <memory>
#include <stddef.h>
#include <string>
#include <vector>

#include "PhredProbabilityTypes.h"
#include "stringFunctions.h"
#include "mathFunctions.h"
#include "probability.h"

#include <armadillo>

namespace GenotypeLikelihoods {
namespace SequencingError {

struct T1stDerivative{
       size_t index;
       double derivative;
};

struct T2ndDerivative{
       size_t index1;
       size_t index2;
       double derivative;
};
//--------------------------------------------------------------
// TCovariateFunction
// Base class for recal covariate functions
//--------------------------------------------------------------
class TFunction {
private:
	size_t _firstParameterIndex;

protected:
	virtual double *_begin() noexcept              = 0;
	virtual double *_end() noexcept                = 0;
	virtual const double *_cbegin() const noexcept = 0;
	virtual const double *_cend() const noexcept   = 0;
	virtual double *_obegin() noexcept             = 0;
	virtual double *_oend() noexcept               = 0;

	void _initializeValues(const std::vector<std::string> &betas);

public:
	TFunction(uint16_t FirstParameterIndex = 0) : _firstParameterIndex(FirstParameterIndex) {}
	virtual ~TFunction() = default;

	// non-virtuals
	void proposeNewParameters(const arma::mat &JxF, uint16_t &index, double lambda) noexcept;
	void rejectProposedParameters() noexcept;
	constexpr size_t firstParameterIndex() const noexcept { return _firstParameterIndex; };

	// virtuals
	virtual uint16_t numParameters() const noexcept               = 0;
	virtual uint16_t numNonZeroFirstDerivatives() const noexcept  = 0;
	virtual uint16_t numNonZeroSecondDerivatives() const noexcept = 0;

	virtual T1stDerivative get1stDerivatives(uint16_t val, size_t i) const noexcept           = 0;
	virtual T2ndDerivative get2ndDerivatives(uint16_t val, size_t i, size_t j) const noexcept = 0;

	// check value range: to ensure that data can be recalibrated
	virtual bool checkOrInitValueRange(const std::vector<uint16_t> &values) = 0;

	// estimation
	virtual double getEtaTerm(uint16_t val) const noexcept                = 0;
	virtual double adjustParametersPostEstimation() noexcept              = 0;
	virtual std::string typeString() const noexcept                       = 0;
	std::string modelString() const;
};

//--------------------------------------------------------------
// TCovariateFunction_intercept
// An intercept term
//--------------------------------------------------------------
class TIntercept : public TFunction {
private:
	double _beta    = 0.;
	double _oldBeta = 0.;

protected:
	double *_begin() noexcept override { return &_beta; }
	double *_end() noexcept override { return &_beta + 1; }
	const double *_cbegin() const noexcept override { return &_beta; }
	const double *_cend() const noexcept override { return &_beta + 1; }
	double *_obegin() noexcept override { return &_oldBeta; }
	double *_oend() noexcept override { return &_oldBeta + 1; }

public:
	static inline const std::string name = "intercept";

	TIntercept(uint16_t FirstParameterIndex, const std::string &beta)
		: TFunction(FirstParameterIndex), _beta(coretools::str::convertStringCheck<double>(beta)){};

	uint16_t numParameters() const noexcept override { return 1; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 1; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	bool checkOrInitValueRange(const std::vector<uint16_t> &) noexcept override { return true; };

	constexpr double intercept() const noexcept { return _beta; }
	constexpr double &intercept() noexcept { return _beta; }

	double getEtaTerm(uint16_t = 0) const noexcept override { return _beta; };
	T1stDerivative get1stDerivatives() const noexcept { return {firstParameterIndex(), 1.}; };
	T1stDerivative get1stDerivatives(uint16_t, size_t) const noexcept override { return {firstParameterIndex(), 1.}; };
	T2ndDerivative get2ndDerivatives(uint16_t, size_t, size_t ) const noexcept override {
		return {firstParameterIndex(), firstParameterIndex(), 0.};
	};
	double adjustParametersPostEstimation() noexcept override { return 0.; }

	virtual std::string typeString() const noexcept override { return name; };
};

//--------------------------------------------------------------
// TCovariateFunction_polynomial
// A polynomial function
//--------------------------------------------------------------
template<size_t O, bool transform = true> class TPolynomial : public TFunction {
	static_assert(O > 0);

private:
	std::array<double, O> _betas{1.};  // betas of the model
	std::array<double, O> _oldBetas{}; // use during estimation

	double _getAsDouble(uint16_t val) const noexcept {
		if constexpr (transform) {
			static const std::array<double, 256> map = []() {
				std::array<double, 256> fs;
				for (size_t v = 0; v < fs.size(); ++v) {
					fs[v] = logit(coretools::Probability(genometools::PhredIntProbability(v)));
				}
				return fs;
			}();
			return map[val];
		} else
			return static_cast<double>(val);
	}

protected:
	double *_begin() noexcept override { return _betas.data(); }
	double *_end() noexcept override { return _betas.data() + O; }
	const double *_cbegin() const noexcept override { return _betas.data(); }
	const double *_cend() const noexcept override { return _betas.data() + O; }
	double *_obegin() noexcept override { return _oldBetas.data(); }
	double *_oend() noexcept override { return _oldBetas.data() + O; }

public:
	static inline const std::string name = "polynomial";

	TPolynomial(uint16_t FirstParameterIndex, const std::vector<std::string> &betas) : TFunction(FirstParameterIndex) {
		_initializeValues(betas);
	}

	uint16_t numParameters() const noexcept override { return O; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return O; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	bool checkOrInitValueRange(const std::vector<uint16_t> &values) noexcept override {
		if constexpr (transform) {
			return std::all_of(values.begin(), values.end(),
							   [](auto v) { return v <= genometools::PhredIntProbability::max().get(); });
		} else
			return true;
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

	T1stDerivative get1stDerivatives(uint16_t val, size_t i) const noexcept override {
		const double v = _getAsDouble(val);
		if constexpr (O == 1) { return {firstParameterIndex(), v}; }
		else return {firstParameterIndex() + i, coretools::uPow(v, i + 1)};
	}

	T2ndDerivative get2ndDerivatives(uint16_t, size_t, size_t) const noexcept override {
		return {firstParameterIndex(), firstParameterIndex(), 0.};
	};

	std::string typeString() const noexcept override {
		using coretools::str::toString;
		return name + '(' + toString(O) + ')';
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

	// tmp storage
	mutable std::vector<TProbitTmpStorage> _tmpStorage;

	void _expandTmpStorage(uint16_t MaxValue) const;

protected:
	double *_begin() noexcept override { return _betas.data(); }
	double *_end() noexcept override { return _betas.data() + _betas.size(); }
	const double *_cbegin() const noexcept override { return _betas.data(); }
	const double *_cend() const noexcept override { return _betas.data() + _betas.size(); }
	double *_obegin() noexcept override { return _oldBetas.data(); }
	double *_oend() noexcept override { return _oldBetas.data() + _oldBetas.size(); }

public:
	static inline const std::string name = "probit";
	TProbit(uint16_t FirstParameterIndex, const std::vector<std::string> &betas);

	uint16_t numParameters() const noexcept override { return 3; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 3; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 6; };

	T1stDerivative get1stDerivatives(uint16_t val, size_t i) const noexcept override;
	T2ndDerivative get2ndDerivatives(uint16_t val, size_t i, size_t j) const noexcept override;

	bool checkOrInitValueRange(const std::vector<uint16_t> &) noexcept override { return true; };

	double getEtaTerm(uint16_t val) const noexcept override;
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

	void _resize(uint16_t MaxValue);
	bool _checkValueRange(uint16_t val) const noexcept { return val < numParameters(); };
	void _adjustValueRanges(const std::vector<uint16_t> &values);

protected:
	double *_begin() noexcept override { return _betas.data(); }
	double *_end() noexcept override { return _betas.data() + _betas.size(); }
	const double *_cbegin() const noexcept override { return _betas.data(); }
	const double *_cend() const noexcept override { return _betas.data() + _betas.size(); }
	double *_obegin() noexcept override { return _oldBetas.data(); }
	double *_oend() noexcept override { return _oldBetas.data() + _betas.size(); }

public:
	static inline const std::string name = "specific";

	TSpecific(uint16_t FirstParameterIndex, const std::vector<std::string> &betas);

	uint16_t numParameters() const noexcept override { return _betas.size(); };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 1; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	T1stDerivative get1stDerivatives(uint16_t val, size_t i) const noexcept override;
	T2ndDerivative get2ndDerivatives(uint16_t, size_t, size_t ) const noexcept override {
		return {firstParameterIndex(), firstParameterIndex(), 0.};
	};

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
	std::vector<double> _betas;            // betas of the model
	std::vector<double> _oldBetas;         // use during estimation
	std::vector<TIndexMapEntry> _indexMap; // maps value to parameter index

	void _resize(size_t NumParameters);
	void _initMapFromVector(const std::vector<uint16_t> &values);
	bool _checkValueRange(uint16_t val) const noexcept { return val < _indexMap.size() && _indexMap[val].used; };
	void _adjustValueRanges(const std::vector<uint16_t> &values);

protected:
	double *_begin() noexcept override { return _betas.data(); }
	double *_end() noexcept override { return _betas.data() + _betas.size(); }
	const double *_cbegin() const noexcept override { return _betas.data(); }
	const double *_cend() const noexcept override { return _betas.data() + _betas.size(); }
	double *_obegin() noexcept override { return _oldBetas.data(); }
	double *_oend() noexcept override { return _oldBetas.data() + _betas.size(); }

public:
	static inline const std::string name = "map";
	TSpecificMap(uint16_t FirstParameterIndex, const std::vector<std::string> &betas);

	uint16_t numParameters() const noexcept override { return _betas.size(); };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 1; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	T1stDerivative get1stDerivatives(uint16_t val, size_t i) const noexcept override;

	T2ndDerivative get2ndDerivatives(uint16_t, size_t, size_t ) const noexcept override {
		return {firstParameterIndex(), firstParameterIndex(), 0.};
	};

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
	virtual std::string typeString() const noexcept override { return name; }
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_ */
