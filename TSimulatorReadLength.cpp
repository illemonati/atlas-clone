/*
 * TSimulatorReadLength.cpp
 *
 *  Created on: Oct 6, 2017
 *      Author: vivian
 */

#include "TSimulatorReadLength.h"


//---------------------------------------------------------
//TSimulatorReadLength
//---------------------------------------------------------
TSimulatorReadLength::TSimulatorReadLength(TRandomGenerator* RandomGenerator, std::string & s){
	randomGenerator = RandomGenerator;

	//is a fixed length
	meanLength = stringToDouble(s);
	if(meanLength < 5 || meanLength > 10000)
		throw "Read length must be between 5 and 10,000!";
	cumulAtMin = 0.0;
};
TSimulatorReadLength::TSimulatorReadLength(TRandomGenerator* RandomGenerator){
	randomGenerator = RandomGenerator;
	meanLength = -1;
	cumulAtMin = 0.0;
};

void TSimulatorReadLength::sample(readLengthContainer & rl){
	rl.fragmentLength = meanLength;
	rl.readLength = meanLength;
};



TSimulatorReadLengthGamma::TSimulatorReadLengthGamma(TRandomGenerator* RandomGenerator, std::string & s):TSimulatorReadLength(RandomGenerator){
	parseFunctionString(s, alpha, beta);
	if(alpha <= 0.0)
		throw "Shape parameter alpha must be > 0.0!";
	if(beta <= 0.0)
		throw "Rate parameter alpha must be > 0.0!";

	calculateAverageLength();
};

TSimulatorReadLengthGamma::TSimulatorReadLengthGamma(TRandomGenerator* RandomGenerator):TSimulatorReadLength(RandomGenerator){
	alpha = -1.0;
	beta = -1.0;
	_min = -1.0;
	_max = -1.0;
}

void TSimulatorReadLengthGamma::parseFunctionString(std::string & s, double & param1, double & param2){
	if(s[0] != '(')
		throw "1 Fail to understand read length function: use format function(var1,var2)[min,max].";
	s.erase(0,1);

	unsigned int pos = s.find(",");
	if(pos == std::string::npos)
		throw "2 Fail to understand read length function: use format function(var1,var2)[min,max].";
	param1 = stringToDouble(s.substr(0,pos));

	s.erase(0,pos+1);

	pos = s.find(")");
	if(pos == std::string::npos)
		throw "3 Fail to understand read length function: use format function(var1,var2)[min,max].";
	param2 = stringToDouble(s.substr(0,pos));
	s.erase(0,pos+1);

	if(s[0] != '[')
		throw "4 Fail to understand read length function: use format function(var1,var2)[min,max].";
	s.erase(0,1);
	pos = s.find(",");
	if(pos == std::string::npos)
		throw "5 Fail to understand read length function: use format function(var1,var2)[min,max].";
	_min = stringToDouble(s.substr(0,pos));
	if(_min <= 0) throw "min read length must be > 0!";
	s.erase(0,pos+1);
	pos = s.find("]");
	if(pos == std::string::npos)
		throw "6 Fail to understand read length function: use format function(var1,var2)[min,max].";
	_max = stringToDouble(s.substr(0,pos));
};

void TSimulatorReadLengthGamma::calculateAverageLength(){
	//get weighted average
	double w;
	double sum = 0.0;
	meanLength = 0.0;
	for(int i=_min; i<_max; ++i){
		w = exp(randomGenerator->gammaLogDensityFunction(i, alpha, beta));
		meanLength += w * (double) i;
		sum += w;
	}
	//add area >= max
	w = 1.0 - randomGenerator->gammaCumulativeDistributionFunction(_max - 0.5, alpha, beta);
	meanLength += w * _max;
	sum += w;
	//normalize
	meanLength /= sum;
	//cumulative at _min
	cumulAtMin = randomGenerator->gammaCumulativeDistributionFunction(_min-0.5, alpha, beta);
}

void TSimulatorReadLengthGamma::sample(readLengthContainer & rl){
	rl.fragmentLength = round(randomGenerator->getGammaRand(alpha, beta));
	while(rl.fragmentLength < _min)
		rl.fragmentLength = round(randomGenerator->getGammaRand(alpha, beta));
	rl.readLength = std::min(rl.fragmentLength, _max);
};

TSimulatorReadLengthGammaMode::TSimulatorReadLengthGammaMode(TRandomGenerator* RandomGenerator, std::string & s):TSimulatorReadLengthGamma(RandomGenerator){
	//here the parameters parsed are mode and variance, so adjust alpha and beta
	parseFunctionString(s, mode, var);
	if(mode <= 0.0)
		throw "Mode of gamma distribution must be > 0.0!";
	if(var <= 0.0)
		throw "Variance of gamma distribution must be > 0.0!";

	beta = (mode + sqrt(mode*mode + 4.0*var)) / (2.0 * var);
	alpha = mode*beta + 1.0;

	if(alpha <= 0.0)
		throw "Provided mode and variance imply a shape parameter alpha <= 0.0!";
	if(beta <= 0.0)
		throw "Provided mode and variance imply a rate parameter beta <= 0.0!";

	calculateAverageLength();
}
