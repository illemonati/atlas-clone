/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Apr 27, 2021
 *      Author: phaentu
 */

#include "TSimulatorQuality.h"

#include <stddef.h>
#include <stdint.h>
#include <algorithm>
#include <memory>

#include "TLog.h"
#include "TRandomGenerator.h"
#include "algorithmsAndVectors.h"
#include "mathFunctions.h"
#include "stringFunctions.h"
#include "strongTypes.h"
#include "weakTypes.h"

namespace Simulations {
using genometools::PhredIntProbability;
using coretools::instances::randomGenerator;
using coretools::str::toString;
using coretools::instances::logfile;

//----------------------------------
// TSimulatorQualityDist
//----------------------------------
TSimulatorQualityDistFixed::TSimulatorQualityDistFixed(std::string &s) {
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

void TSimulatorQualityDistFixed::printDetails(const std::string &Name) const {
	logfile().list(Name + ": fixed quality of " + toString(_max));
}

//------------------------------------------------
// TSimulatorQualityDistBinned
//------------------------------------------------
TSimulatorQualityDistBinned::TSimulatorQualityDistBinned(std::string &s){
	const auto pos1 = s.find("(");
	if (pos1 == 0) {
		s.erase(0, 1);
		const auto pos2 = s.find(')');
		if (pos2 == std::string::npos || pos2 != s.size() - 1)
			throw "Failed to understand binned distribution '" + s + "'! Use binned(quality_1,quality_2,..,quality_n).";
		s.erase(pos2, 1);
		coretools::str::fillContainerFromString(s, _qualBins, ',', false, true, false);
	} else {
		throw "Failed to understand binned distribution '" + s + "'! Use binned(quality_1,quality_2,..,quality_n).";
	}
}

PhredIntProbability TSimulatorQualityDistBinned::sample() const noexcept {
	return _qualBins[randomGenerator().sample(_qualBins.size())];
}

void TSimulatorQualityDistBinned::printDetails(const std::string &Name) const {
	logfile().list(Name + ": uniformly distributed among the values " + coretools::str::concatenateString(_qualBins, ", "));
}

//------------------------------------------------
// TSimulatorQualityDistFreq
//------------------------------------------------
TSimulatorQualityDistFreq::TSimulatorQualityDistFreq(std::string &s){
	const auto pos1 = s.find("(");
	if (pos1 == 0) {
		s.erase(0, 1);
		const auto pos2 = s.find(')');
		if (pos2 == std::string::npos || pos2 != s.size() - 1) {
			throw "Failed to understand frequency distribution '" + s +
				"'! Use frequency(quality_1:frequency_1,quality_2:frequency_2, ... ,quality_n:frequency_n).";
		}
		s.erase(pos2, 1);
		std::vector<std::string> tmp;
		coretools::str::fillContainerFromString(s, tmp, ',');

		// now parse each bin
		_qualBins.resize(tmp.size());
		_frequencies.resize(tmp.size());

		for (size_t i = 0; i < tmp.size(); ++i) {
			const auto pos3 = tmp[i].find(':');
			if (pos3 == std::string::npos) {
				throw "Failed to understand frequency distribution '" + s +
					"'! Use frequency(quality_1:frequency_1,quality_2:frequency_2, ... ,quality_n:frequency_n).";
			}
			_qualBins[i]    = tmp[i].substr(0, pos3);
			_frequencies[i] = tmp[i].substr(pos3 + 1);
		}

		// fill cumulative
		coretools::fillCumulative(_frequencies, _cumulativeFrequencies);

	} else {
		throw "Failed to understand frequency distribution '" + s +
			"'! Use frequency(quality_1:frequency_1,quality_2:frequency_2, ... ,quality_n:frequency_n).";
	}
}

PhredIntProbability TSimulatorQualityDistFreq::sample() const noexcept {
	return _qualBins[randomGenerator().pickOne(_cumulativeFrequencies)];
}

void TSimulatorQualityDistFreq::printDetails(const std::string &Name) const {
	logfile().list(Name + ": frequency bins " +
	       coretools::str::concatenateString(coretools::str::paste(_qualBins, _frequencies, ":"), ", "));
}

//------------------------------------------------
// TSimulatorQualityDistNormal
//------------------------------------------------
TSimulatorQualityDistNormal::TSimulatorQualityDistNormal(std::string &s) {
	parseFunctionString(s);
	fillDensities();
}

TSimulatorQualityDistNormal::TSimulatorQualityDistNormal(double mean, double sd, int min, int max)
	: _mean(mean), _sd(sd), _min(min), _max(max) {
	fillDensities();
}

void TSimulatorQualityDistNormal::parseFunctionString(std::string &s) {
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
	_min = coretools::str::convertString<PhredIntProbability>(s.substr(0, pos));
	s.erase(0, pos + 1);
	pos = s.find("]");
	if (pos == std::string::npos)
		throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
	_max = coretools::str::convertString<PhredIntProbability>(s.substr(0, pos));
	if (_max < _min) throw "Fail to understand distribution '" + orig + "': max must be >= min!";
}

void TSimulatorQualityDistNormal::fillDensities() {
	using namespace coretools::TNormalDistr;
	// fill densities
	const auto size = _max.get() + 1 - _min.get();
	_densities.resize(size);
	_cumulDensities.resize(size);

	double nextDens = cumulativeDistrFunction(_min.get() - 0.5, _mean, _sd * _sd);
	double prevDens;
	double sum = 0;
	for (int i = 0; i < size; ++i) {
		prevDens      = nextDens;
		nextDens      = cumulativeDistrFunction(_min.get() + i + 0.5, _mean, _sd * _sd);
		_densities[i] = nextDens - prevDens;
		sum += _densities[i];
	}

	// normalize
	for (int i = 0; i < size; ++i) _densities[i] /= sum;

	// fill cumulative
	_cumulDensities[0] = _densities[0];
	for (int i = 1; i < size; ++i) _cumulDensities[i] = _cumulDensities[i - 1] + _densities[i];
	_cumulDensities[size - 1] = 1.0;
}

PhredIntProbability TSimulatorQualityDistNormal::sample() const noexcept {
	return PhredIntProbability(randomGenerator().pickOne(_cumulDensities) + _min.get());
}

void TSimulatorQualityDistNormal::printDetails(const std::string &Name) const {
	logfile().list(Name + ": Normally distributed quality scores with mean=", _mean, " and sd=", _sd,
		      ", truncated to [", _min, ",", _max, "].");
}

} // namespace Simulations
