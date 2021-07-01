/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Apr 27, 2021
 *      Author: phaentu
 */

#include "TSimulatorQuality.h"

namespace Simulations{

//----------------------------------
//TSimulatorQualityDist
//----------------------------------
TSimulatorQualityDist::TSimulatorQualityDist(std::string & s, TRandomGenerator* RandomGenerator){
	size_t pos = s.find("(");
	if(pos == std::string::npos)
		_max = coretools::str::convertStringCheck<uint8_t>(s);
	else if(pos == 0){
		pos = s.find(')');
		if(pos == std::string::npos || pos != s.size() - 1)
			throw "Failed to understand fixed distribution '" + s + "'!";
		_max = coretools::str::convertStringCheck<uint8_t>(s.substr(1,pos - 1));
	} else {
		throw "Failed to understand fixed distribution '" + s + "'!";
	}

	_min = 0;
	_maxPlusOne = -1;
	_mean = -1.0;
	_sd = -1.0;
	_randomGenerator = RandomGenerator;
};

TSimulatorQualityDist::TSimulatorQualityDist(TRandomGenerator* RandomGenerator){
	_max = 30;
	_min = 0;
	_maxPlusOne = -1;
	_mean = -1.0;
	_sd = -1.0;
	_randomGenerator = RandomGenerator;
};

void TSimulatorQualityDist::sample(std::vector<PhredIntProbability> & phredInt) const{
	for(auto& q : phredInt){
		q = sample();
	}
};

std::string TSimulatorQualityDist::_details() const{
	return "fixed quality of " + toString(_max);
};

void TSimulatorQualityDist::printDetails(coretools::TLog* logfile, const std::string & Name) const{
	logfile->list(Name + ": " + _details());
};

//------------------------------------------------
// TSimulatorQualityDistBinned
//------------------------------------------------
TSimulatorQualityDistBinned::TSimulatorQualityDistBinned(std::string & s, TRandomGenerator* RandomGenerator):TSimulatorQualityDist(RandomGenerator){
	size_t pos = s.find("(");
	if(pos == 0){
		s.erase(0,1);
		pos = s.find(')');
		if(pos == std::string::npos || pos != s.size() - 1)
			throw "Failed to understand binned distribution '" + s + "'! Use binned(quality_1,quality_2,..,quality_n).";
		s.erase(pos,1);
		coretools::str::fillContainerFromString(s, _qualBins, ',', false, true, false);
	} else {
		throw "Failed to understand binned distribution '" + s + "'! Use binned(quality_1,quality_2,..,quality_n).";
	}
};

PhredIntProbability TSimulatorQualityDistBinned::sample() const{
	size_t tmp = _randomGenerator->sample(_qualBins.size());
	return _qualBins[tmp];
};

std::string TSimulatorQualityDistBinned::_details() const{
	std::string tmpS = coretools::str::concatenateString(_qualBins, ", ");
	return "uniformly distributed among the values " + tmpS;
};

//------------------------------------------------
// TSimulatorQualityDistFreq
//------------------------------------------------
TSimulatorQualityDistFreq::TSimulatorQualityDistFreq(std::string & s, TRandomGenerator* RandomGenerator):TSimulatorQualityDist(RandomGenerator){
	size_t pos = s.find('(');
	if(pos == 0){
		s.erase(0,1);
		pos = s.find(')');
		if(pos == std::string::npos || pos != s.size() - 1){
			throw "Failed to understand frequency distribution '" + s + "'! Use frequency(quality_1:frequency_1,quality_2:frequency_2, ... ,quality_n:frequency_n).";
		}
		s.erase(pos,1);
		std::vector<std::string> tmp;
		coretools::str::fillContainerFromString(s, tmp, ',');

		//now parse each bin
		_qualBins.resize(tmp.size());
		_frequencies.resize(tmp.size());

		for(size_t i = 0; i < tmp.size(); ++i){
			pos = tmp[i].find(':');
			if(pos == std::string::npos){
				throw "Failed to understand frequency distribution '" + s + "'! Use frequency(quality_1:frequency_1,quality_2:frequency_2, ... ,quality_n:frequency_n).";
			}
			_qualBins[i] = tmp[i].substr(0, pos);
			_frequencies[i] = tmp[i].substr(pos+1);
		}

		//fill cumulative
		coretools::fillCumulative(_frequencies, _cumulativeFrequencies);

	} else {
		throw "Failed to understand frequency distribution '" + s + "'! Use frequency(quality_1:frequency_1,quality_2:frequency_2, ... ,quality_n:frequency_n).";
	}
};

PhredIntProbability TSimulatorQualityDistFreq::sample() const{
	return _qualBins[_randomGenerator->pickOne(_cumulativeFrequencies)];
};

std::string TSimulatorQualityDistFreq::_details() const{
	return "frequency bins " + coretools::str::concatenateString(coretools::str::paste(_qualBins, _frequencies, ":"), ", ");
};

//------------------------------------------------
// TSimulatorQualityDistNormal
//------------------------------------------------
TSimulatorQualityDistNormal::TSimulatorQualityDistNormal(std::string & s, TRandomGenerator* RandomGenerator):TSimulatorQualityDist(RandomGenerator){
	parseFunctionString(s);
	_maxPlusOne = _max.get() + 1;

	fillDensities();
};

TSimulatorQualityDistNormal::TSimulatorQualityDistNormal(double & mean, double & sd, int min, int max, TRandomGenerator* RandomGenerator):TSimulatorQualityDist(RandomGenerator){
	_mean = mean;
	_sd = sd;
	_min = min;
	_max = max;
	_maxPlusOne = _max.get() + 1;

	fillDensities();
};

void TSimulatorQualityDistNormal::parseFunctionString(std::string & s){
	std::string orig = s;

	if(s[0] != '(')
		throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
	s.erase(0,1);

	unsigned int pos = s.find(",");
	if(pos == std::string::npos)
		throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
	_mean = coretools::str::convertString<double>(s.substr(0,pos));
	if(_mean < 0)
		throw "Fail to understand distribution '" + orig + "': mean must be > 0.";
	s.erase(0,pos+1);

	pos = s.find(")");
	if(pos == std::string::npos)
		throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
	_sd = coretools::str::convertString<double>(s.substr(0,pos));
	if(_sd < 0)
			throw "Fail to understand distribution '" + orig + "': sd must be > 0.";
	s.erase(0,pos+1);

	if(s[0] != '[')
		throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
	s.erase(0,1);
	pos = s.find(",");
	if(pos == std::string::npos)
		throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
	_min = coretools::str::convertString<PhredIntProbability>(s.substr(0,pos));
	s.erase(0,pos+1);
	pos = s.find("]");
	if(pos == std::string::npos)
		throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
	_max = coretools::str::convertString<PhredIntProbability>(s.substr(0,pos));
	if(_max < _min)
			throw "Fail to understand distribution '" + orig + "': max must be >= min!";
};

void TSimulatorQualityDistNormal::fillDensities(){
	//fill densities
	_size = _maxPlusOne.get() - _min.get();
	_densities.resize(_size);
	_cumulDensities.resize(_size);

	double nextDens = coretools::TNormalDistr::cumulativeDistrFunction(_min.get() - 0.5, _mean, _sd*_sd);
	double prevDens;
	double sum = 0;
	for(int i=0; i<_size; ++i){
		prevDens = nextDens;
		nextDens = coretools::TNormalDistr::cumulativeDistrFunction(_min.get() + i + 0.5, _mean, _sd*_sd);
		_densities[i] =  nextDens - prevDens;
		sum += _densities[i];
	}

	//normalize
	for(int i=0; i<_size; ++i)
		_densities[i] /= sum;

	//fill cumulative
	_cumulDensities[0] = _densities[0];
	for(int i=1; i<_size; ++i)
		_cumulDensities[i] = _cumulDensities[i-1] + _densities[i];
	_cumulDensities[_size-1] = 1.0;
}

PhredIntProbability TSimulatorQualityDistNormal::sample() const {
	return PhredIntProbability( _randomGenerator->pickOne(_cumulDensities) + _min.get() );
};

std::string TSimulatorQualityDistNormal::_details() const{
	return coretools::str::toString("Normally distributed quality scores with mean=", _mean, " and sd=", _sd, ", truncated to [", _min, ",", _max, "].");
};

}; //end namesapce

