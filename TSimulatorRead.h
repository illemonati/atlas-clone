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
private:
	TLog* logfile;
	TRandomGenerator* randomGenerator;
	char* toBase;

	//PMD
	TPMD* pmdObject;
	bool pmdInitialized = false;

	//qual params
	double mQ;
	double sdQ;
	double meanQual, sdQual;
	int maxQual;
	bool qualToErroTableInitialized = false;
	double* qualToErroTable;

	//final read params
	std::string queryBases = "";
	std::string bamQualities = "";

	//helper functions
	double dePhred(double x);

	//general functions
	void setQualityDistribution(double mean, double sd, int maxQ);
	int sampleQuality();
	void initializeQualToErrorTable();

	//virtual functions
	virtual int returnQual(int qual, int pos, BaseContext baseContext, int maxQual);

protected:
	void applyPMD(short & base, long & posInRead, readLengthContainer & rl);

public:
	TSimulatorRead(TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase);
	virtual ~TSimulatorRead(){
		if(qualToErroTableInitialized)
			delete[] qualToErroTable;
	};

	virtual void simulate(short* posAddress, readLengthContainer & rl, TGenotypeMap & genoMap);

	//getters and setters
	std::string getQueryBases();
	std::string getBamQualities();

};


//-------------------------------
// recal transformation
//-------------------------------

class TSimulatorReadRecal:public TSimulatorRead{
private:
	std::vector<double> betas;
	double* qualTermForTransformation;
	double* posTermForTransformation;

	//private functions
	int transformQuality(int & qual, int & pos, int & context);

//	TSimulatorRecalTransform(std::vector<double> & Betas, TSimulatorReadLength& readLengthDist);

	//derived functions
	int returnQual(int & qual, int & pos, int & context);

public:
	TSimulatorReadRecal(std::vector<double> Betas, int & maxReadLen, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase);
	~TSimulatorReadRecal(){};
};


//-------------------------------
// BQSR transformation
//-------------------------------

class TSimulatorBQSRTransform:public TSimulatorRead{
	float alpha;
	float beta;
	std::string readGroupName;

public:
	TSimulatorBQSRTransform(int & maxReadLen, TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
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
