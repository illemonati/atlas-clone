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
	TRandomGenerator* _randomGenerator;
	uint32_t _meanLength;
	double _cumulAtMin;

	double* _positionProbs; //normalized (1 - cumulDensity)

public:
	TReadLengthDistribution(std::string & s, TRandomGenerator* RandomGenerator);
	TReadLengthDistribution(TRandomGenerator* RandomGenerator);
	virtual ~TReadLengthDistribution();

	double operator[](const uint32_t position) const{ return _positionProbs[position]; };

	virtual TReadLength sample();
	virtual uint32_t max(){return _meanLength;};
	virtual double mean(){return _meanLength;};
	virtual double probAcceptance(){return 1.0 - _cumulAtMin;};
	virtual void printDetails(TLog* logfile);

};

class TSimulatorReadLengthGamma:public TReadLengthDistribution{
protected:
	double _meanLength;
	double _alpha, _beta;
	uint32_t _min, _maxPlusOne;
	bool _initialized;

	double* _gammaDensity;
	double* _gammaCumulDensity;

	void parseFunctionString(std::string & s, double & param1, double & param2);
	void initiate(TLog* logfile);

public:
	TSimulatorReadLengthGamma(std::string & s, TRandomGenerator* RandomGenerator, TLog* Logfile);
	TSimulatorReadLengthGamma(TRandomGenerator* RandomGenerator);
	virtual ~TSimulatorReadLengthGamma(){
		if(_initialized){
			delete[] _gammaDensity;
			delete[] _gammaCumulDensity;
			delete[] _positionProbs;
		}
	};
	TReadLength sample();
	virtual uint32_t max(){return _maxPlusOne - 1;};
	virtual double mean(){return _meanLength;};
	virtual void printDetails(TLog* logfile);
};

class TSimulatorReadLengthGammaMode:public TSimulatorReadLengthGamma{
protected:
	double _mode, _var;

public:
	TSimulatorReadLengthGammaMode(std::string & s, TRandomGenerator* RandomGenerator, TLog* Logfile);
	~TSimulatorReadLengthGammaMode(){};
	void printDetails(TLog* logfile);
};

}; //end namespace

#endif /* TSIMULATORREADLENGTH_H_ */
