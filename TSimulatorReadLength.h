/*
 * TSimulatorReadLength.h
 *
 *  Created on: Oct 6, 2017
 *      Author: vivian
 */

#ifndef TSIMULATORREADLENGTH_H_
#define TSIMULATORREADLENGTH_H_

#include "commonutilities/stringFunctions.h"
#include "commonutilities/TRandomGenerator.h"
#include "commonutilities/TLog.h"

//---------------------------------------------------------
//TSimulatorReadLength
//---------------------------------------------------------
class TSimulatorReadLength{
protected:
	TRandomGenerator* randomGenerator;
	int meanLength;
	double cumulAtMin;

public:
	TSimulatorReadLength(std::string & s, TRandomGenerator* RandomGenerator);
	TSimulatorReadLength(TRandomGenerator* RandomGenerator);
	double* positionProbs; //normalized (1 - cumulDensity)
	virtual ~TSimulatorReadLength(){};

	virtual void sample(int & readLength, int & fragmentLength);
	virtual int max(){return meanLength;};
	virtual double mean(){return meanLength;};
	virtual double probAcceptance(){return 1.0 - cumulAtMin;};
	virtual void printDetails(TLog* logfile);

};

class TSimulatorReadLengthGamma:public TSimulatorReadLength{
protected:
	double meanLength;
	double alpha, beta;
	int _min, _max;
	bool initialized;

	double* gammaDensity;
	double* gammaCumulDensity;

	void parseFunctionString(std::string & s, double & param1, double & param2);
	void initiate();

public:
	TSimulatorReadLengthGamma(std::string & s, TRandomGenerator* RandomGenerator);
	TSimulatorReadLengthGamma(TRandomGenerator* RandomGenerator);
	virtual ~TSimulatorReadLengthGamma(){
		if(initialized){
			delete[] gammaDensity;
			delete[] gammaCumulDensity;
			delete[] positionProbs;
		}
	};
	void sample(int & readLength, int & fragmentLength);
	virtual int max(){return _max;};
	virtual double mean(){return meanLength;};
	virtual void printDetails(TLog* logfile);
};

class TSimulatorReadLengthGammaMode:public TSimulatorReadLengthGamma{
protected:
	double mode, var;

public:
	TSimulatorReadLengthGammaMode(std::string & s, TRandomGenerator* RandomGenerator);
	~TSimulatorReadLengthGammaMode(){};
	void printDetails(TLog* logfile);
};



#endif /* TSIMULATORREADLENGTH_H_ */
