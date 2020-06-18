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
#include "TLog.h"

namespace Simulations{

//---------------------------------------------------------
//ReadLength
//---------------------------------------------------------
struct TReadLength{
	uint16_t read;
	uint16_t fragment;
	TReadLength(const uint16_t Read, const uint16_t Fragment){
		read = Read;
		fragment = Fragment;
	};

	uint16_t diff(){
		return fragment - read;
	};
};

//---------------------------------------------------------
//TReadLengthDistribution
//---------------------------------------------------------
class TReadLengthDistribution{
protected:
	TRandomGenerator* randomGenerator;
	int meanLength;
	double cumulAtMin;

public:
	double* positionProbs; //normalized (1 - cumulDensity)

	TReadLengthDistribution(std::string & s, TRandomGenerator* RandomGenerator);
	TReadLengthDistribution(TRandomGenerator* RandomGenerator);
	virtual ~TReadLengthDistribution();

	virtual TReadLength sample();
	virtual int max(){return meanLength;};
	virtual double mean(){return meanLength;};
	virtual double probAcceptance(){return 1.0 - cumulAtMin;};
	virtual void printDetails(TLog* logfile);

};

class TSimulatorReadLengthGamma:public TReadLengthDistribution{
protected:
	double meanLength;
	double alpha, beta;
	uint16_t _min, _maxPlusOne;
	bool initialized;

	double* gammaDensity;
	double* gammaCumulDensity;

	void parseFunctionString(std::string & s, double & param1, double & param2);
	void initiate(TLog* logfile);

public:
	TSimulatorReadLengthGamma(std::string & s, TRandomGenerator* RandomGenerator, TLog* Logfile);
	TSimulatorReadLengthGamma(TRandomGenerator* RandomGenerator);
	virtual ~TSimulatorReadLengthGamma(){
		if(initialized){
			delete[] gammaDensity;
			delete[] gammaCumulDensity;
			delete[] positionProbs;
		}
	};
	TReadLength sample();
	virtual int max(){return _maxPlusOne - 1;};
	virtual double mean(){return meanLength;};
	virtual void printDetails(TLog* logfile);
};

class TSimulatorReadLengthGammaMode:public TSimulatorReadLengthGamma{
protected:
	double mode, var;

public:
	TSimulatorReadLengthGammaMode(std::string & s, TRandomGenerator* RandomGenerator, TLog* Logfile);
	~TSimulatorReadLengthGammaMode(){};
	void printDetails(TLog* logfile);
};

}; //end namespace

#endif /* TSIMULATORREADLENGTH_H_ */
