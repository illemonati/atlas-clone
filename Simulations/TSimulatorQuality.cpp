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
TSimulatorQualityDist::TSimulatorQualityDist(std::string & s){
	size_t pos = s.find("(");
	if(pos == std::string::npos)
		_max = convertStringCheck<int>(s);
	else if(pos == 0){
		pos = s.find(')');
		if(pos == std::string::npos || pos != s.size() - 1)
			throw "Failed to understand fixed distribution '" + s + "'!";
		_max = convertStringCheck<int>(s.substr(1,pos - 1));
	} else
		throw "Failed to understand fixed distribution '" + s + "'!";

	_min = -1;
	_maxPlusOne = -1;
	_mean = -1.0;
	_sd = -1.0;
};

TSimulatorQualityDist::TSimulatorQualityDist(){
	_max = 30;
	_min = -1;
	_maxPlusOne = -1;
	_mean = -1.0;
	_sd = -1.0;
}

void TSimulatorQualityDist::sample(int* qualities, const int & len){
	for(int i = 0; i < len; ++i)
		qualities[i] = _max;
};

void TSimulatorQualityDist::printDetails(TLog* logfile, const std::string & Name) const{
	logfile->list(Name + ": fixed quality of " + toString(_max));
};

//------------------------------------------------
// TSimulatorQualityDistBinned
//------------------------------------------------
TSimulatorQualityDistBinned::TSimulatorQualityDistBinned(std::string & s, TRandomGenerator* RandomGenerator):TSimulatorQualityDist(){
	size_t pos = s.find("(");
	if(pos == 0){
		s.erase(0,1);
		pos = s.find(')');
		if(pos == std::string::npos || pos != s.size() - 1)
			throw "Failed to understand binned distribution '" + s + "'! Use binned(quality_1,quality_2,..,quality_n).";
		s.erase(pos,1);
		fillVectorFromString(s, _qualBins, ',');
	} else {
		throw "Failed to understand binned distribution '" + s + "'! Use binned(quality_1,quality_2,..,quality_n).";
	}

	_randomGenerator = RandomGenerator;
};

void TSimulatorQualityDistBinned::sample(int* qualities, const int & len){
	for(int i=0; i<len; ++i){
		qualities[i] = _qualBins[_randomGenerator->sample(_qualBins.size())];
	}
};

void TSimulatorQualityDistBinned::printDetails(TLog* logfile, const std::string & Name) const{
	std::string tmpS = concatenateString(_qualBins, ", ");
	logfile->list(Name + ": uniformly distributed among the values " + tmpS);
};

//------------------------------------------------
// TSimulatorQualityDistFreq
//------------------------------------------------
TSimulatorQualityDistFreq::TSimulatorQualityDistFreq(std::string & s, TRandomGenerator* RandomGenerator):TSimulatorQualityDist(){
	size_t pos = s.find('(');
	if(pos == 0){
		s.erase(0,1);
		pos = s.find(')');
		if(pos == std::string::npos || pos != s.size() - 1){
			throw "Failed to understand frequency distribution '" + s + "'! Use frequency(quality_1:frequency_1,quality_2:frequency_2, ... ,quality_n:frequency_n).";
		}
		s.erase(pos,1);
		std::vector<std::string> tmp;
		fillVectorFromString(s, tmp, ',');

		//now parse each bin
		_qualBins.resize(tmp.size());
		_frequencies.resize(tmp.size());

		for(size_t i = 0; i < tmp.size(); ++i){
			pos = tmp[i].find(':');
			if(pos == std::string::npos){
				throw "Failed to understand frequency distribution '" + s + "'! Use frequency(quality_1:frequency_1,quality_2:frequency_2, ... ,quality_n:frequency_n).";
			}
			convertString(tmp[i].substr(0, pos), _qualBins[i]);
			convertString(tmp[i].substr(pos+1), _frequencies[i]);
		}

		//fill cumulative
		fillCumulative(_frequencies, _cumulativeFrequencies);

	} else {
		throw "Failed to understand frequency distribution '" + s + "'! Use frequency(quality_1:frequency_1,quality_2:frequency_2, ... ,quality_n:frequency_n).";
	}

	_randomGenerator = RandomGenerator;
};

void TSimulatorQualityDistFreq::sample(int* qualities, const int & len){
	for(int i=0; i<len; ++i){
		qualities[i] = _qualBins[_randomGenerator->pickOne(_cumulativeFrequencies)];
	}
};

void TSimulatorQualityDistFreq::printDetails(TLog* logfile, const std::string & Name) const{
	std::string tmpS = concatenateString(_qualBins, ",");
	logfile->list(Name + ": frequency bins " + concatenateString(paste(_qualBins, _frequencies, ":"), ", "));
};

//------------------------------------------------
// TSimulatorQualityDistNormal
//------------------------------------------------
TSimulatorQualityDistNormal::TSimulatorQualityDistNormal(std::string & s, TRandomGenerator* RandomGenerator):TSimulatorQualityDist(){
	_randomGenerator = RandomGenerator;

	parseFunctionString(s);
	_maxPlusOne = _max + 1;

	fillDensities();
};

TSimulatorQualityDistNormal::TSimulatorQualityDistNormal(double & mean, double & sd, int min, int max, TRandomGenerator* RandomGenerator):TSimulatorQualityDist(){
	_randomGenerator = RandomGenerator;

	_mean = mean;
	_sd = sd;
	_min = min;
	_max = max;
	_maxPlusOne = _max + 1;

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
	_mean = convertString<double>(s.substr(0,pos));
	if(_mean < 0)
		throw "Fail to understand distribution '" + orig + "': mean must be > 0.";
	s.erase(0,pos+1);

	pos = s.find(")");
	if(pos == std::string::npos)
		throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
	_sd = convertString<double>(s.substr(0,pos));
	if(_sd < 0)
			throw "Fail to understand distribution '" + orig + "': sd must be > 0.";
	s.erase(0,pos+1);

	if(s[0] != '[')
		throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
	s.erase(0,1);
	pos = s.find(",");
	if(pos == std::string::npos)
		throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
	_min = convertString<double>(s.substr(0,pos));
	if(_min < 0)
		throw "Fail to understand function '" + orig + "': min must be >= 0!";
	s.erase(0,pos+1);
	pos = s.find("]");
	if(pos == std::string::npos)
		throw "Fail to understand distribution '" + orig + "': use format normal(mean,sd)[min,max].";
	_max = convertString<double>(s.substr(0,pos));
	if(_max < _min)
			throw "Fail to understand distribution '" + orig + "': max must be >= min!";
};

void TSimulatorQualityDistNormal::fillDensities(){
	//fill densities
	_size = _maxPlusOne - _min;
	_densities.resize(_size);
	_cumulDensities.resize(_size);

	double nextDens = TNormalDistr::cumulativeDistrFunction(_min-0.5, _mean, _sd*_sd);
	double prevDens;
	double sum = 0;
	for(int i=0; i<_size; ++i){
		prevDens = nextDens;
		nextDens = TNormalDistr::cumulativeDistrFunction(_min + i + 0.5, _mean, _sd*_sd);
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

int TSimulatorQualityDistNormal::sample(){
	return _randomGenerator->pickOne(_cumulDensities) + _min;
};

void TSimulatorQualityDistNormal::sample(int* qualities, const int & len){
	for(uint16_t tmpInt=0; tmpInt<len; ++tmpInt){
		qualities[tmpInt] = _randomGenerator->pickOne(_cumulDensities) + _min;
	}
};

void TSimulatorQualityDistNormal::printDetails(TLog* logfile){
	logfile->list("Normally distributed quality scores with mean=" + toString(_mean) + " and sd=" + toString(_sd) + ", truncated to [" + toString(_min) + "," + toString(_max) + "].");
};


}; //end namesapce

