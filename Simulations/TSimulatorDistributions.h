/*
 * TSimulatorQuality.h
 *
 *  Created on: Apr 27, 2021
 *      Author: phaentu
 */

#ifndef SIMULATIONS_TSIMULATORDISTRIBUTIONS_H_
#define SIMULATIONS_TSIMULATORDISTRIBUTIONS_H_

#include <string>
#include <vector>
#include "TLog.h"
#include "TRandomGenerator.h"
#include "probability.h"
#include "algorithms.h"

namespace Simulations {

using coretools::instances::logfile;
using coretools::instances::randomGenerator;

// Pure abstract base class
template <typename valueType>
class TSimulatorDistribution {
protected:
	TSimulatorDistribution() = default;
public:
	virtual ~TSimulatorDistribution() = default;
	virtual void sample(std::vector<valueType> & vectorToFill) const noexcept final {
		for (auto &q : vectorToFill){
			q = sample();
		}
	}
	virtual valueType sample() const noexcept                                = 0;
	virtual void printDetails(const std::string &Name) const = 0;
};

//-------------------------------
// TSimulatorDistributionFixed
// Used for quality and mapping quality
//-------------------------------
// Class of a fixed value
template <typename valueType>
class TSimulatorDistributionFixed : public TSimulatorDistribution<valueType> {
private:
	valueType _max{30};
public:
	TSimulatorDistributionFixed(std::string &s){
		const auto pos1 = s.find("(");
		if (pos1 == std::string::npos)
			_max = coretools::str::convertStringCheck<uint8_t>(s);
		else if (pos1 == 0) {
			const auto pos2 = s.find(')');
			if (pos2 != s.size() - 1)
				throw "Failed to understand fixed distribution '" + s + "'!";
			_max = coretools::str::convertStringCheck<uint8_t>(s.substr(1, pos2 - 1));
		} else {
			throw "Failed to understand fixed distribution '" + s + "'!";
		}
	}

	valueType sample() const noexcept override { return _max; };
	void printDetails(const std::string &Name) const override {
		logfile().list(Name, ": fixed value of ", _max);
	}
};

//------------------------------------------------
// TSimulatorDistributionBinned
// Class of a binned distribution
//------------------------------------------------
template <typename valueType>
class TSimulatorDistributionBinned : public TSimulatorDistribution<valueType> {
private:
	std::vector<valueType> _qualBins;
public:
	TSimulatorDistributionBinned(std::string &s){
		const auto pos1 = s.find("(");
		if (pos1 == 0) {
			s.erase(0, 1);
			const auto pos2 = s.find(')');
			if (pos2 == std::string::npos || pos2 != s.size() - 1)
				throw "Failed to understand binned distribution '" + s + "'! Use binned(quality_1,quality_2,..,quality_n).";
			s.erase(pos2, 1);
			coretools::str::fillContainerFromString(s, _qualBins, ',', true, true, false);
		} else {
			throw "Failed to understand binned distribution '" + s + "'! Use binned(quality_1,quality_2,..,quality_n).";
		}
	}

	valueType sample() const noexcept override{
		return _qualBins[randomGenerator().sample(_qualBins.size())];
	}
	void printDetails(const std::string &Name) const override{
		logfile().list(Name + ": uniformly distributed among the values " + coretools::str::concatenateString(_qualBins, ", "));
	}
};

//------------------------------------------------
// TSimulatorDistributionFreq
// Class of a binned distribution with frequencies
//------------------------------------------------
template <typename valueType>
class TSimulatorDistributionFreq : public TSimulatorDistribution<valueType> {
private:
	std::vector<valueType> _valueBins;
	std::vector<coretools::Probability> _frequencies;
	std::vector<coretools::Probability> _cumulativeFrequencies;

public:
	TSimulatorDistributionFreq(std::string &s){
	const auto pos1 = s.find("(");
		if (pos1 == 0) {
			s.erase(0, 1);
			const auto pos2 = s.find(')');
			if (pos2 == std::string::npos || pos2 != s.size() - 1) {
				throw "Failed to understand frequency distribution '" + s +
					"'! Use frequency(value_1:frequency_1,value_2:frequency_2, ... ,value_n:frequency_n).";
			}
			s.erase(pos2, 1);
			std::vector<std::string> tmp;
			coretools::str::fillContainerFromString(s, tmp, ',');

			// now parse each bin
			_valueBins.resize(tmp.size());
			_frequencies.resize(tmp.size());

			for (size_t i = 0; i < tmp.size(); ++i) {
				const auto pos3 = tmp[i].find(':');
				if (pos3 == std::string::npos) {
					throw "Failed to understand frequency distribution '" + s +
						"'! Use frequency(value_1:frequency_1,value_2:frequency_2, ... ,value_n:frequency_n).";
				}
				_valueBins[i]    = tmp[i].substr(0, pos3);
				_frequencies[i] = tmp[i].substr(pos3 + 1);
			}

			// fill cumulative
			coretools::fillCumulative(_frequencies, _cumulativeFrequencies);

		} else {
			throw "Failed to understand frequency distribution '" + s +
				"'! Use frequency(value_1:frequency_1,value_2:frequency_2, ... ,value_n:frequency_n).";
		}
	}

	valueType sample() const noexcept override{
		return _valueBins[randomGenerator().pickOne(_cumulativeFrequencies)];
	}
	void printDetails(const std::string &Name) const override{
		logfile().list(Name + ": frequency bins " +
			       coretools::str::concatenateString(coretools::str::paste(_valueBins, _frequencies, ":"), ", "));
	}
};

//------------------------------------------------
// TSimulatorDistributionDiscretized
// Base class for all distributions that discretize continuous distributions
//------------------------------------------------
template <typename valueType>
class TSimulatorDistributionDiscretized : public TSimulatorDistribution<valueType> {
protected:
	// densities
	std::vector<double> _cumulDensities;
	valueType _truncationMin;
	valueType _truncationMax;

public:
	TSimulatorDistributionDiscretized() = default;
	TSimulatorDistributionDiscretized(valueType TruncationMin, valueType TruncationMax): _truncationMin(TruncationMin), _truncationMax(TruncationMax){};

	valueType sample() const noexcept override{
		return randomGenerator().pickOne(_cumulDensities) + _truncationMin.get();
	}
};

//------------------------------------------------
// TSimulatorQualityDistNormal
// Class of a normal distribution
//------------------------------------------------

template <typename valueType>
class TSimulatorDistributionNormal : public TSimulatorDistributionDiscretized<valueType> {
private:
	// densities
	std::vector<double> _cumulDensities;
	double _mean = -1.;
	double _sd   = -1.;
	using TSimulatorDistributionDiscretized<valueType>::_truncationMin;
	using TSimulatorDistributionDiscretized<valueType>::_truncationMax;

	void _parseFunctionString(std::string &s){
		std::string orig = s;

		if (s[0] != '(') throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
		s.erase(0, 1);

		auto pos = s.find(",");
		if (pos == std::string::npos)
			throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
		_mean = coretools::str::convertString<double>(s.substr(0, pos));
		if (_mean < 0) throw "Fail to understand distribution '" + orig + "': mean must be > 0.";
		s.erase(0, pos + 1);

		pos = s.find(")");
		if (pos == std::string::npos)
			throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
		_sd = coretools::str::convertString<double>(s.substr(0, pos));
		if (_sd < 0) throw "Fail to understand distribution '" + orig + "': sd must be > 0.";
		s.erase(0, pos + 1);

		if (s[0] != '[') throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
		s.erase(0, 1);
		pos = s.find(",");
		if (pos == std::string::npos)
			throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
		_truncationMin = coretools::str::convertString<valueType>(s.substr(0, pos));
		s.erase(0, pos + 1);
		pos = s.find("]");
		if (pos == std::string::npos)
			throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
		_truncationMax = coretools::str::convertString<valueType>(s.substr(0, pos));
		if (_truncationMax < _truncationMin) throw "Fail to understand distribution '" + orig + "': max must be >= min!";
	}

	void _fillDensities(){
		using namespace coretools::TNormalDistr;
		// fill densities
		//const auto size = _truncationMax.get() + 1 - _truncationMax.get();
		const size_t size = _truncationMax - _truncationMax + 1;
		std::vector<double> densities(size);
		_cumulDensities.resize(size);

		double nextDens = cumulativeDistrFunction(_truncationMin.get() - 0.5, _mean, _sd * _sd);
		double prevDens;
		for (size_t i = 0; i < size; ++i) {
			prevDens      = nextDens;
			nextDens      = cumulativeDistrFunction(_truncationMin.get() + i + 0.5, _mean, _sd * _sd);
			densities[i] = nextDens - prevDens;
		}

		// fill cumulative
		coretools::fillCumulative(densities, _cumulDensities);
	}

public:
	TSimulatorDistributionNormal(std::string &s){
		_parseFunctionString(s);
		_fillDensities();
	}

	TSimulatorDistributionNormal(double Mean, double SD, valueType TruncationMin, valueType TruncationMax):
		TSimulatorDistributionDiscretized<valueType>(TruncationMin, TruncationMax), _mean(Mean), _sd(SD) {
		_fillDensities();
	}

	void printDetails(const std::string &Name) const override{
		logfile().list(Name + ": Normally distributed with mean=", _mean, " and sd=", _sd,
				      ", truncated to [", _truncationMin, ",", _truncationMax, "].");
	}
};

} // namespace Simulations

#endif /* SIMULATIONS_TSIMULATORDISTRIBUTIONS_H_ */
