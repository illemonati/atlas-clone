
/*
 * TPMDFunction.h
 *
 *  Created on: Feb 23, 2023
 *      Author: andreas
 */

#ifndef TPMDFUNCTION_H_
#define TPMDFUNCTION_H_

#include <armadillo>
#include <map>
#include <string>

#include "TPMDTables.h"
#include "coretools/Types/probability.h"

namespace GenotypeLikelihoods {
using TPMDEstimationParameters = std::map<std::string, double>;

class TPMDFunction {
public:
	TPMDFunction()          = default;
	virtual ~TPMDFunction() = default;

	virtual bool hasDamage() const noexcept = 0;

	virtual void parseEstimationParameters(TPMDEstimationParameters &EstimationParameters) = 0;
	virtual void learn(const TPMDTable &Table, const genometools::Base &from, const genometools::Base &to,
					   const TPMDEstimationParameters &EstimationParameters)               = 0;
	virtual std::string string() const noexcept                                            = 0;
	virtual coretools::Probability prob(uint16_t pos) const noexcept                       = 0;
};

class TPMDFunctionNoPMD final : public TPMDFunction {
public:
	static inline const std::string name    = "none";
	static inline const std::string example = name;
	TPMDFunctionNoPMD(std::string_view string);

	bool hasDamage() const noexcept override { return false; }
	std::string string() const noexcept override { return name + "[]"; }

	void parseEstimationParameters(TPMDEstimationParameters &) override{};
	void learn(const TPMDTable &, const genometools::Base &, const genometools::Base &,
			   const TPMDEstimationParameters &) override{};

	coretools::Probability prob(uint16_t) const noexcept override { return 0.0; }
};

class TPMDFunctionExponential final : public TPMDFunction {
private:
	uint16_t _lastPosition;
	double _a, _b, _c;
	std::vector<coretools::Probability> _values;

	void _initialEstimatesOLS(const countVec &pmdCounts, const countVec &pmdSums, std::array<double, 3> &Parameters);
	void _fillFAndJacobian(arma::vec &F, arma::mat &J, const countVec &pmdCounts, const countVec &pmdSums,
						   const std::array<double, 3> &Parameters);
	void _estimateWithNewtonRaphson(const countVec &pmdCounts, const countVec &pmdSums,
									std::array<double, 3> &Parameters, uint32_t numNRIterations, double epsilon);
	double _calcLL(const countVec &pmdCounts, const countVec &pmdSums, const std::array<double, 3> &Parameters);
	void _fillPMDProbabilities();

public:
	static inline const std::string name    = "Exponential";
	static inline const std::string example = name + "[lastPosition,a,b,c]";
	static inline const std::string epsilon = "PMDExponentialEpsilon";
	static inline const std::string numNR   = "PMDExponentialNumNR";
	TPMDFunctionExponential(std::string_view string);

	bool hasDamage() const noexcept override { return _lastPosition > 0; }
	std::string string() const noexcept override {
		return name + "[" + coretools::str::toString(_lastPosition) + ',' +
			   coretools::str::concatenateString(std::vector{_a, _b, _c}, ",") + "]";
	}

	void parseEstimationParameters(TPMDEstimationParameters &EstimationParameters) override;
	void learn(const TPMDTable &Table, const genometools::Base &from, const genometools::Base &to,
			   const TPMDEstimationParameters &EstimationParameters) override;

	coretools::Probability prob(uint16_t pos) const noexcept override;
};

class TPMDFunctionEmpiric final : public TPMDFunction {
private:
	std::vector<coretools::Probability> _values;

public:
	static inline const std::string name    = "Empiric";
	static inline const std::string example = name + "[p1,p2,...]";
	TPMDFunctionEmpiric(std::string_view string);

	bool hasDamage() const noexcept override { return _values.size() + _values.back().get() != 1.0; }
	std::string string() const noexcept override {
		return name + "[" + coretools::str::concatenateString(_values, ",") + "]";
	}

	void parseEstimationParameters(TPMDEstimationParameters &) override{};
	void learn(const TPMDTable &Table, const genometools::Base &from, const genometools::Base &to,
			   const TPMDEstimationParameters &EstimationParameters) override;

	coretools::Probability prob(uint16_t pos) const noexcept override;
};

TPMDFunction *makeFunction(std::string_view pmdString);
} // namespace GenotypeLikelihoods

#endif
