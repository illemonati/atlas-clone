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
TSimulatorReadLength::TSimulatorReadLength(std::string & s, TRandomGenerator* RandomGenerator){
	randomGenerator = RandomGenerator;

	//expect string (x) -> remov ( and )!
	s.erase(0, 1);
	s.erase(s.length()-1, 1);
	meanLength = stringToInt(s);
	if(meanLength < 5 || meanLength > 10000)
		throw "Read length must be between 5 and 10,000!";

	cumulAtMin = 0.0;

	positionProbs = new double[meanLength];
	for(int i=0; i<meanLength; ++i){
		positionProbs[i] = 1.0 / meanLength;
	}
};

TSimulatorReadLength::TSimulatorReadLength(TRandomGenerator* RandomGenerator){
	randomGenerator = RandomGenerator;
	meanLength = -1;
	cumulAtMin = 0.0;

	positionProbs = NULL;
};

TSimulatorReadLength::~TSimulatorReadLength(){
	if(meanLength > 0)
		delete[] positionProbs;
}

void TSimulatorReadLength::sample(int & readLength, int & fragmentLength){
	fragmentLength = meanLength;
	readLength = meanLength;
};

void TSimulatorReadLength::printDetails(TLog* logfile){
	logfile->list("Reads of fixed length " + toString(meanLength) + ".");
};

//--------------------------------------------------
//TSimulatorReadLengthGamma
//--------------------------------------------------
TSimulatorReadLengthGamma::TSimulatorReadLengthGamma(std::string & s, TRandomGenerator* RandomGenerator):TSimulatorReadLength(RandomGenerator){
	std::cout << "in gamma cstr with s" << std::endl;
	parseFunctionString(s, alpha, beta);
	if(alpha <= 0.0)
		throw "Shape parameter alpha must be > 0.0!";
	if(beta <= 0.0)
		throw "Rate parameter alpha must be > 0.0!";

	initiate();
};

TSimulatorReadLengthGamma::TSimulatorReadLengthGamma(TRandomGenerator* RandomGenerator):TSimulatorReadLength(RandomGenerator){
	std::cout << "in gamma cstr without s" << std::endl;

	alpha = -1.0;
	beta = -1.0;
	_min = -1.0;
	_max = -1.0;
	meanLength = -1.0;

	gammaDensity = NULL;
	gammaCumulDensity = NULL;
	positionProbs = NULL;

	initialized = false;
}

void TSimulatorReadLengthGamma::parseFunctionString(std::string & s, double & param1, double & param2){
	std::string orig = s;

	if(s[0] != '(')
		throw "Fail to understand function '" + orig + "': use format function(var1,var2)[min,max].";
	s.erase(0,1);

	unsigned int pos = s.find(",");
	if(pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format function(var1,var2)[min,max].";
	param1 = stringToDouble(s.substr(0,pos));

	s.erase(0,pos+1);

	pos = s.find(")");
	if(pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format function(var1,var2)[min,max].";
	param2 = stringToDouble(s.substr(0,pos));
	s.erase(0,pos+1);

	if(s[0] != '[')
		throw "Fail to understand function '" + orig + "': use format function(var1,var2)[min,max].";
	s.erase(0,1);
	pos = s.find(",");
	if(pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format function(var1,var2)[min,max].";
	_min = stringToDouble(s.substr(0,pos));
	if(_min <= 0)
		throw "Fail to understand function '" + orig + "': min read length must be > 0!";
	s.erase(0,pos+1);
	pos = s.find("]");
	if(pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format function(var1,var2)[min,max].";
	_max = stringToDouble(s.substr(0,pos));
	if(_max <= _min)
			throw "Fail to understand function '" + orig + "': max must be > min!";
};

void TSimulatorReadLengthGamma::initiate(){
	//prepare storage
	gammaDensity = new double[_max];
	gammaCumulDensity = new double[_max];
	positionProbs = new double[_max];
	initialized = true;

	//1) calc density and get weighted average
	//first set all bins < _min to zero
	double totalArea = 0.0;
	for(int i=0; i < _min; ++i)	gammaDensity[i] = 0;

	//then calculate densities for all bins <_max
	double gammaDensity_i;
	for(int i=_min; i<(_max-1); ++i){
		gammaDensity[i] = exp(randomGenerator->gammaLogDensityFunction(i, alpha, beta));
		totalArea += gammaDensity[i];
	}

	//add area >= max and 0's for < min to table
	gammaDensity[_max - 1] = (1.0 - randomGenerator->gammaCumulativeDistributionFunction(_max - 0.5, alpha, beta));
	totalArea += gammaDensity[_max - 1];

	//normalize densities (needed because truncated at _min)
	//also calc mean read length
	meanLength = 0.0;
	for(int i=_min; i<_max; ++i){
		gammaDensity[i] /= totalArea;
		meanLength += i * gammaDensity[i];
	}

	//2) make table for cumulative gamma distribution
	gammaCumulDensity[0] = gammaDensity[0];
	for(int i=1; i < _max; ++i)	gammaCumulDensity[i] = gammaCumulDensity[i-1] + gammaDensity[i];

	//3) distribution of position probabilities (=normalized 1 - cumul)
	positionProbs[0] = 1.0; //position 1 is always present in read
	double sum = positionProbs[0];
	for(int i=1; i < _max; ++i){
		positionProbs[i] = 1.0 -  gammaCumulDensity[i-1];
		sum += positionProbs[i];
	}

	//normalize
	for(int i=0; i < _max; ++i)
		positionProbs[i] /= sum;
}

void TSimulatorReadLengthGamma::sample(int & readLength, int & fragmentLength){
	fragmentLength = round(randomGenerator->getGammaRand(alpha, beta));
	while(fragmentLength < _min)
		fragmentLength = round(randomGenerator->getGammaRand(alpha, beta));
	readLength = std::min(fragmentLength, _max);
}

void TSimulatorReadLengthGamma::printDetails(TLog* logfile){
	logfile->list("Gamma distributed fragment length with alpha=" + toString(alpha) + " and beta=" + toString(beta) + " of at least " + toString(_min) + ".");
	logfile->list("Fragments  > " + toString(_max) + " will result in reads of length " + toString(_max) + ".");
	if(probAcceptance() < 0.9)
		logfile->warning("The chosen distribution will only result in " + toString(probAcceptance()) + " of draws being accepted.");
};


//--------------------------------------------------
//TSimulatorReadLengthGammaMode
//--------------------------------------------------
TSimulatorReadLengthGammaMode::TSimulatorReadLengthGammaMode(std::string & s, TRandomGenerator* RandomGenerator):TSimulatorReadLengthGamma(RandomGenerator){
	//here the parameters parsed are mode and variance, so adjust alpha and beta
	std::cout << "in cstr gamma mode" << std::endl;
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

void TSimulatorReadLengthGammaMode::printDetails(TLog* logfile){
	logfile->list("Gamma distributed fragment length with mode=" + toString(mode) + " and variance=" + toString(var) + " of at least " + toString(_min) + ".");
	logfile->list("Fragments  > " + toString(_max) + " will result in reads of length " + toString(_max) + ".");
	if(probAcceptance() < 0.9)
		logfile->warning("The chosen distribution will only result in " + toString(probAcceptance()) + " of draws being accepted.");
};
