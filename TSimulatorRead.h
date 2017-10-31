/*
 * TSimulatorQuality.h
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#ifndef TSIMULATORREAD_H_
#define TSIMULATORREAD_H_

#include "TGenotypeMap.h"
#include "TSimulatorReadLength.h"
#include "TPostMortemDamage.h"


class TSimulatorRead{
public:
	TSimulatorReadLength* readLengthDist;
	TPMD* pmdObject;
	TParameters params;
	TLog* logfile;
	TRandomGenerator* randomGenerator;

	double mQ;
	double sdQ;
	double meanQual, sdQual;
	int maxQual;
	bool pmdInitialized = false;
	bool qualToErroTableInitialized = false;
	double* qualToErroTable;



	TSimulatorRead(TSimulatorReadLength* ReadLengthDist, TPMD* PmdObject);
	int sampleQuality();
	void setQualityTransformation(std::vector<double> Betas);
	void initializeQualToErrorTable();

	virtual ~TSimulatorRead(){
		if(qualToErroTableInitialized)
			delete[] qualToErroTable;
	};
	virtual int returnQual(int qual, int pos, BaseContext baseContext, int maxQual);

};

class TSimulatorReadRecal:public TSimulatorRead{
public:
	double* beta;
	double* qualTermForTransformation;
	double* posTermForTransformation;

	TSimulatorReadRecal(std::vector<double> Betas, TSimulatorReadLength* ReadLengthDist, TPMD* PmdObject);
	int transformQuality(int & qual, int pos, int context);

//	TSimulatorRecalTransform(std::vector<double> & Betas, TSimulatorReadLength& readLengthDist);
	int returnQual(int qual, int pos, BaseContext baseContext, int maxQual);

	~TSimulatorReadRecal(){};

};

class TSimulatorBQSRTransform:public TSimulatorRead{
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


#endif /* TSIMULATORREAD_H_ */
