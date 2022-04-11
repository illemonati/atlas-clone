/*
 * TSimulatorReadLength.cpp
 *
 *  Created on: Oct 6, 2017
 *      Author: vivian
 */

#include "TSimulatorReadLength.h"

#include <math.h>
#include <stddef.h>
#include <algorithm>

#include "TLog.h"
#include "TRandomGenerator.h"
#include "commonWeakTypes.h"
#include "mathFunctions.h"
#include "stringFunctions.h"
#include "weakTypes.h"

namespace Simulations {

using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using coretools::str::convertString;

//---------------------------------------------------------
// TReadLengthDistribution
//---------------------------------------------------------
TReadLengthDistribution::TReadLengthDistribution(std::string &s) {

	// expect string (x) -> remov ( and )!
	s.erase(0, 1);
	s.erase(s.length() - 1, 1);
	_meanLength = convertString<int>(s);
	if (_meanLength < 5 || _meanLength > 10000) throw "Read length must be between 5 and 10,000!";

	_positionProbs.resize(_meanLength, 1./_meanLength);
}

void TReadLengthDistribution::printDetails() {
	logfile().list("Reads of fixed length ", _meanLength, ".");
}

//--------------------------------------------------
// TSimulatorReadLengthGamma
//--------------------------------------------------
TSimulatorReadLengthGamma::TSimulatorReadLengthGamma(std::string &s) {
	parseFunctionString(s, _alpha, _beta);
	if (_alpha <= 0.0) throw "Shape parameter alpha must be > 0.0!";
	if (_beta <= 0.0) throw "Rate parameter alpha must be > 0.0!";
	initiate();
}

void TSimulatorReadLengthGamma::parseFunctionString(std::string &s, double &param1, double &param2) {
	std::string orig = s;

	if (s[0] != '(') throw "Fail to understand function '" + orig + "': use format function(var1,var2)[min,max].";
	s.erase(0, 1);

	auto pos = s.find(",");
	if (pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format function(var1,var2)[min,max].";
	param1 = convertString<double>(s.substr(0, pos));

	s.erase(0, pos + 1);

	pos = s.find(")");
	if (pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format function(var1,var2)[min,max].";
	param2 = convertString<double>(s.substr(0, pos));
	s.erase(0, pos + 1);

	if (s[0] != '[') throw "Fail to understand function '" + orig + "': use format function(var1,var2)[min,max].";
	s.erase(0, 1);
	pos = s.find(",");
	if (pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format function(var1,var2)[min,max].";
	_min = convertString<double>(s.substr(0, pos));
	if (_min <= 0) throw "Fail to understand function '" + orig + "': min read length must be > 0!";
	s.erase(0, pos + 1);
	pos = s.find("]");
	if (pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format function(var1,var2)[min,max].";
	_maxPlusOne = convertString<double>(s.substr(0, pos)) + 1;
	if (_maxPlusOne < _min) throw "Fail to understand function '" + orig + "': max must be > min!";
}

void TSimulatorReadLengthGamma::initiate() {
	using coretools::TGammaDistr;
	// prepare storage
	_gammaDensity.resize(_maxPlusOne, 0.);
	_gammaCumulDensity.resize(_maxPlusOne);
	_positionProbs.resize(_maxPlusOne);

	// 1) calc density and get weighted average
	// first set all bins < _min to zero
	double totalArea = 0.0;

	// then calculate densities for all bins <_max
	for (size_t i = _min; i < (_maxPlusOne - 1); ++i) {
		_gammaDensity[i] = TGammaDistr::density(i, _alpha, _beta);
		totalArea += _gammaDensity[i];
	}

	// add area >= max
	_gammaDensity[_maxPlusOne - 1] =
		1.0 - TGammaDistr::cumulativeDistrFunction(_maxPlusOne - 0.5, _alpha, _beta);
	totalArea += _gammaDensity[_maxPlusOne - 1];

	// normalize densities (needed because truncated at _min)
	// also calc mean read length
	_meanLength = 0;
	for (uint32_t i = _min; i < _maxPlusOne; ++i) {
		_gammaDensity[i] /= totalArea;
		_meanLength += i * _gammaDensity[i];
	}

	// 2) make table for cumulative gamma distribution
	_gammaCumulDensity[0] = _gammaDensity[0];
	for (uint32_t i = 1; i < _maxPlusOne; ++i) _gammaCumulDensity[i] = _gammaCumulDensity[i - 1] + _gammaDensity[i];

	// 3) distribution of position probabilities (=normalized 1 - cumul)
	_positionProbs[0] = 1.0; // position 1 is always present in read
	double sum        = _positionProbs[0];
	for (uint32_t i = 1; i < _maxPlusOne; ++i) {
		_positionProbs[i] = 1.0 - _gammaCumulDensity[i - 1];
		sum += _positionProbs[i];
	}

	// normalize
	for (auto & p: _positionProbs) p/=sum;

	if (_gammaCumulDensity[_min] > 0.5)
		logfile().warning("This readLength distribution results in ", _gammaCumulDensity[_min] * 100,
				 "% discarded fragments because they are smaller than the minimum read length! Choose "
				 "different parameters to reduce run time.");
}

TReadLength TSimulatorReadLengthGamma::sample() const noexcept {
	uint16_t fragmentLength = round(randomGenerator().getGammaRand(_alpha, _beta));
	while (fragmentLength < _min) fragmentLength = round(randomGenerator().getGammaRand(_alpha, _beta));

	if (fragmentLength < _maxPlusOne - 1) {
		return TReadLength(fragmentLength, fragmentLength);
	} 
	return TReadLength(_maxPlusOne - 1, fragmentLength);
}

void TSimulatorReadLengthGamma::printDetails() {
	logfile().list("Gamma distributed fragment length with alpha=", _alpha, " and beta=", _beta, " of at least ", _min,
		      ".");
	logfile().list("Fragments  > ", _maxPlusOne, " will result in reads of length ", _maxPlusOne, ".");
	if (probAcceptance() < 0.9)
		logfile().warning("The chosen distribution will only result in ", probAcceptance(), " of draws being accepted.");
}

//--------------------------------------------------
// TSimulatorReadLengthGammaMode
//--------------------------------------------------
TSimulatorReadLengthGammaMode::TSimulatorReadLengthGammaMode(std::string &s) : TSimulatorReadLengthGamma() {
	// here the parameters parsed are mode and variance, so adjust alpha and beta
	parseFunctionString(s, _mode, _var);
	if (_mode <= 0.0) throw "Mode of gamma distribution must be > 0.0!";
	if (_var <= 0.0) throw "Variance of gamma distribution must be > 0.0!";

	_beta  = (_mode + sqrt(_mode * _mode + 4.0 * _var)) / (2.0 * _var);
	_alpha = _mode * _beta + 1.0;

	if (_alpha <= 0.0) throw "Provided mode and variance imply a shape parameter alpha <= 0.0!";
	if (_beta <= 0.0) throw "Provided mode and variance imply a rate parameter beta <= 0.0!";

	initiate();
}

void TSimulatorReadLengthGammaMode::printDetails() {
	logfile().list("Gamma distributed fragment length with mode=", _mode, " and variance=", _var, " of at least ", _min,
		      ".");
	logfile().list("Fragments  > ", _maxPlusOne, " will result in reads of length ", _maxPlusOne, ".");
	if (probAcceptance() < 0.9)
		logfile().warning("The chosen distribution will only result in ", probAcceptance(), " of draws being accepted.");
};

}; // namespace Simulations
