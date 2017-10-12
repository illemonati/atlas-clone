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
	virtual ~TSimulatorQuality(){};
	virtual int returnQual(int qual, int pos, BaseContext baseContext, int maxQual);

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
	int returnQual(int qual, int pos, BaseContext baseContext, int maxQual);

	~TSimulatorRecalTransform(){};

};

class TSimulatorBQSRTransform:public TSimulatorQuality{
	float alpha;
	float beta;
	std::string readGroupName;

public:
	TSimulatorBQSRTransform(std::string qualTransform, TSimulatorReadLength* readLengthDist);
	virtual ~TSimulatorBQSRTransform(){};
	int returnQual(int qual, int pos, BaseContext baseContext, int maxQual);

};

class TSimulatorBQSRPositionTransform:public TSimulatorBQSRTransform{
private:
	float revIntercept;
	float transformationSlope;
	float calculateSum(float & curSlope, float & curIntercept, float & newSum);
	void findTransformationSlope();
public:
	TSimulatorBQSRPositionTransform(float positionTransform, std::string qualTransform, TSimulatorReadLength* readLengthDist);
	~TSimulatorBQSRPositionTransform(){};
};


#endif /* TSIMULATORQUALITY_H_ */
