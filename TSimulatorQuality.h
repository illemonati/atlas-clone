/*
 * TSimulatorQuality.h
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#ifndef TSIMULATORQUALITY_H_
#define TSIMULATORQUALITY_H_

#include "TGenotypeMap.h"
#include "TSimulatorReadLength.h"


class TSimulatorQuality{
public:
	TSimulatorReadLength* readLengthDist;

	TSimulatorQuality(TSimulatorReadLength* ReadLengthDist);
	TSimulatorQuality(){};
	virtual ~TSimulatorQuality(){};
	virtual int returnQual(int qual, int pos, BaseContext baseContext);

};

class TSimulatorRecalTransform:public TSimulatorQuality{
public:
	double* beta;
	double* qualTermForTransformation;
	double* posTermForTransformation;

	TSimulatorRecalTransform(std::vector<double> Betas, TSimulatorReadLength* ReadLengthDist);
	void setQualityTransformation(std::vector<double> Betas);
	int transformQuality(int & qual, int pos, int context);

//	TSimulatorRecalTransform(std::vector<double> & Betas, TSimulatorReadLength& readLengthDist);
	int returnQual(int qual, int pos, BaseContext baseContext);

	~TSimulatorRecalTransform(){};

};



#endif /* TSIMULATORQUALITY_H_ */
