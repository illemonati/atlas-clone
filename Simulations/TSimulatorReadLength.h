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
//TSimulatorReadLength
//---------------------------------------------------------
class TSimulatorReadLength{
protected:
	TRandomGenerator* randomGenerator;
	int meanLength;
	double cumulAtMin;

public:
	double* positionProbs; //normalized (1 - cumulDensity)

	TSimulatorReadLength(std::string & s, TRandomGenerator* RandomGenerator);
	TSimulatorReadLength(TRandomGenerator* RandomGenerator);
	virtual ~TSimulatorReadLength();

	virtual void sample(uint16_t & readLength, uint16_t & fragmentLength);
	virtual int max(){return meanLength;};
	virtual double mean(){return meanLength;};
	virtual double probAcceptance(){return 1.0 - cumulAtMin;};
	virtual void printDetails(TLog* logfile);

};

class TSimulatorReadLengthGamma:public TSimulatorReadLength{
protected:
	double meanLength;
	double alpha, beta;
	int _min, _maxPlusOne;
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
	void sample(uint16_t & readLength, uint16_t & fragmentLength);
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
