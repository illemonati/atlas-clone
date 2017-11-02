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


//toBase[0] = 'A'; toBase[1] = 'C'; toBase[2] = 'G'; toBase[3] = 'T';



class TSimulatorRead{
private:
	//PMD
	TPMD* pmdObject;
	bool pmdInitialized = false;

	TLog* logfile;
	TRandomGenerator* randomGenerator;
	TGenotypeMap genoMap;

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
	char toBase[4] = {'A', 'C', 'G', 'T'};
	double dePhred(double x);

public:
	TSimulatorRead(TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	virtual ~TSimulatorRead(){
		if(qualToErroTableInitialized)
			delete[] qualToErroTable;
	};

	//getters and setters
	std::string getQueryBases();
	std::string getBamQualities();

	//general functions
//	void setPMD(TPMD* PmdObject);
	void setQualityDistribution(double mean, double sd, int maxQ);
	int sampleQuality();

	//virtual functions
	void simulate(short* posAddress, readLengthContainer & rl);
	void setQualityTransformation(std::vector<double> Betas);
	void initializeQualToErrorTable();
	virtual int returnQual(int qual, int pos, BaseContext baseContext, int maxQual);

	//helper tools

};

class TSimulatorReadRecal:public TSimulatorRead{
public:
	double* beta;
	double* qualTermForTransformation;
	double* posTermForTransformation;

	TSimulatorReadRecal(std::vector<double> Betas, int & maxReadLen, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	int transformQuality(int & qual, int & pos, int & context);

//	TSimulatorRecalTransform(std::vector<double> & Betas, TSimulatorReadLength& readLengthDist);
	int returnQual(int & qual, int & pos, int & context);

	~TSimulatorReadRecal(){};

};

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
