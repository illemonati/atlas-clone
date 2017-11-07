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
	meanLength = stringToInt(s);
	if(meanLength < 5 || meanLength > 10000)
		throw "Read length must be between 5 and 10,000!";

	gammaDensity = new float[meanLength + 1];
	gammaCumulDensity = new float[meanLength + 1];
	for(int i=0; i<(meanLength); ++i){
		gammaDensity[i] = 0.0;
		gammaCumulDensity[i] = 0.0;
	}
	gammaDensity[meanLength - 1] = 1.0;
	gammaCumulDensity[meanLength] = 1.0;
	cumulAtMin = 0.0;
};

TSimulatorReadLength::TSimulatorReadLength(TRandomGenerator* RandomGenerator){
	randomGenerator = RandomGenerator;
	meanLength = -1;
	cumulAtMin = 0.0;
	gammaDensity = NULL;
	gammaCumulDensity = NULL;
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

	initiate();
};

TSimulatorReadLengthGamma::TSimulatorReadLengthGamma(TRandomGenerator* RandomGenerator):TSimulatorReadLength(RandomGenerator){
	alpha = -1.0;
	beta = -1.0;
	_min = -1.0;
	_max = -1.0;
	meanLength = -1.0;
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

void TSimulatorReadLengthGamma::initiate(){
	cumulAtMin = randomGenerator->gammaCumulativeDistributionFunction(_min-0.5, alpha, beta);
	double totalArea = 1 - cumulAtMin;

	gammaDensity = new float[_max];
	gammaCumulDensity = new float[_max];

	//get weighted average
	double w;
	double sum = 0.0;
	meanLength = 0.0;
	for(int i=_min; i<_max; ++i){
		double gammaDensity_i = exp(randomGenerator->gammaLogDensityFunction(i, alpha, beta));
		w = gammaDensity_i;
		meanLength += w * (double) i;
		sum += w;

		//fill gammaDensity table at the same time
		gammaDensity[i] = gammaDensity_i;
	}

	//add area >= max and 0's for < min to table
	for(int i=0; i < _min; ++i)	gammaDensity[i] = 0;
	gammaDensity[_max - 1] = 1.0 - randomGenerator->gammaCumulativeDistributionFunction(_max - 0.5, alpha, beta) / totalArea;

	//add area >= max to average calculation
	w = 1.0 - randomGenerator->gammaCumulativeDistributionFunction(_max - 0.5, alpha, beta);
	meanLength += w * _max;
	sum += w;
	meanLength /= sum; 	//normalize

	//make table for cumulative gamma distribution
	for(int i=1; i < _max; ++i){
		gammaCumulDensity[i] = gammaDensity[i-1] + gammaDensity[i];
	}
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

	initiate();
}
