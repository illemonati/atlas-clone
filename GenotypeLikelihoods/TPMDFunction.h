
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
#include <vector>

#include "coretools/Types/probability.h"

namespace GenotypeLikelihoods {

class TPMDFunction {
public:
	TPMDFunction()          = default;
	virtual ~TPMDFunction() = default;

	virtual bool hasDamage() const noexcept                                                    = 0;
	virtual void learn(const std::vector<double> &from_to, const std::vector<double> &to_from) = 0;
	virtual std::string string() const noexcept                                                = 0;
	virtual coretools::Probability prob(uint16_t pos) const noexcept                           = 0;
};

class TPMDFunctionNoPMD final : public TPMDFunction {
public:
	static inline const std::string name    = "none";
	static inline const std::string example = name;
	TPMDFunctionNoPMD(std::string_view string);

	bool hasDamage() const noexcept override { return false; }
	std::string string() const noexcept override { return name + "[]"; }
	void learn(const std::vector<double> &, const std::vector<double> &) override {}
	coretools::Probability prob(uint16_t) const noexcept override { return 0.0; }
};

class TPMDFunctionExponential final : public TPMDFunction {
private:
	double _a, _b, _c;
	std::vector<coretools::Probability> _values;
	void _fillPMDProbabilities(size_t N);

public:
	static inline const std::string name    = "Exponential";
	static inline const std::string example = name + "[lastPosition,a,b,c]";
	TPMDFunctionExponential(std::string_view string);

	bool hasDamage() const noexcept override { return _values.size() > 0; }
	std::string string() const noexcept override {
		return name + "[" + coretools::str::toString(_values.size()) + ',' +
			   coretools::str::concatenateString(std::vector{_a, _b, _c}, ",") + "]";
	}

	void learn(const std::vector<double> &from_to, const std::vector<double> &to_from) override;
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

	void learn(const std::vector<double> &from_to, const std::vector<double> &to_from) override;
	coretools::Probability prob(uint16_t pos) const noexcept override;
};

TPMDFunction *makeFunction(std::string_view pmdString);
} // namespace GenotypeLikelihoods

#endif
