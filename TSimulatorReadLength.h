/*
 * TSimulatorReadLength.h
 *
 *  Created on: Oct 6, 2017
 *      Author: vivian
 */

#ifndef TSIMULATORREADLENGTH_H_
#define TSIMULATORREADLENGTH_H_

#include "stringFunctions.h"
#include "TRandomGenerator.h"


//---------------------------------------------------------
//TSimulatorReadLength
//---------------------------------------------------------
struct readLengthContainer{
	int fragmentLength;
	int readLength;
};

class TSimulatorReadLength{
protected:
	TRandomGenerator* randomGenerator;
	int meanLength;
	double cumulAtMin;
	float* gammaDensity;
	float* gammaCumulDensity;

public:
	TSimulatorReadLength(TRandomGenerator* RandomGenerator, std::string & s);
	TSimulatorReadLength(TRandomGenerator* RandomGenerator);
	virtual ~TSimulatorReadLength(){
		delete[] gammaDensity;
		delete[] gammaCumulDensity;
	};

	virtual void sample(readLengthContainer & rl);
	virtual int max(){return meanLength;};
	virtual double mean(){return meanLength;};
	virtual double probAcceptance(){return 1.0 - cumulAtMin;};
	virtual std::string getFunctionString(){ return "Will simulate reads of fixed length " + toString(meanLength) + ".";};
};

class TSimulatorReadLengthGamma:public TSimulatorReadLength{
protected:
	double alpha, beta;
	int _min, _max;

	void parseFunctionString(std::string & s, double & param1, double & param2);
	void initiate();

public:
	TSimulatorReadLengthGamma(TRandomGenerator* RandomGenerator, std::string & s);
	TSimulatorReadLengthGamma(TRandomGenerator* RandomGenerator);
	virtual ~TSimulatorReadLengthGamma(){};
	void sample(readLengthContainer & rl);
	virtual int max(){return _max;};
	virtual std::string getFunctionString(){ return "Will simulate reads of gamma distributed length with alpha=" + toString(alpha) + " and beta=" + toString(beta) + ".";};
};

class TSimulatorReadLengthGammaMode:public TSimulatorReadLengthGamma{
protected:
	double mode, var;

public:
	TSimulatorReadLengthGammaMode(TRandomGenerator* RandomGenerator, std::string & s);
	virtual ~TSimulatorReadLengthGammaMode(){};
	std::string getFunctionString(){ return "Will simulate reads of gamma distributed length with mode=" + toString(mode) + " and variance=" + toString(var) + ".";};
};



#endif /* TSIMULATORREADLENGTH_H_ */
