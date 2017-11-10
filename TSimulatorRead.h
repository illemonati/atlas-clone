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

protected:
	TLog* logfile;
	TRandomGenerator* randomGenerator;
	char* toBase;

	//PMD
	TPMD* pmdObject;
	bool pmdInitialized = false;

	//final read params
	std::string queryBases = "";
	std::string bamQualities = "";

	//helper functions
//	double dePhredAscii(int x);
	double dePhred(int x);
	int phred(double x);

	//virtual functions
	virtual int returnBamQual(int qual, int pos, BaseContext baseContext, int maxQual);

	//qual params
	int minQual = 0;
	bool dePhredTableInitialized = false;
	double* dePhredTable;
	double meanQual, sdQual;
	int maxQual, maxQualPlusOne;
	bool qualToErroTableInitialized = false;
	double* qualToErroTable;

	//general functions
	void setTrueQualityDistribution(double mean, double sd, int maxQ);
	virtual int sampleTrueQuality();
	void initializeDePhredTable();
	void initializeQualToErrorTable();
	void applyPMD(short & base, long & posInRead, readLengthContainer & rl);

public:
	TSimulatorRead(TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase);
	virtual ~TSimulatorRead(){
		if(qualToErroTableInitialized)
			delete[] qualToErroTable;
		if(dePhredTableInitialized)
			delete[] dePhredTable;
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

	//derived functions
//	int returnQual(int & qual, int & pos, int & context);

public:
	TSimulatorReadRecal(std::vector<double> Betas, int & maxReadLen, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase);
	virtual ~TSimulatorReadRecal(){};
};

//-------------------------------
// BQSR transformation pos
//-------------------------------

class TSimulatorReadBQSR:public TSimulatorRead{
private:

protected:
	TSimulatorReadLength* readLengthDist;
	int maxPos;

	//quality params
	int phi1;
	double phi2;
	bool trueQualToFakeQualTableInitialized = false;
	double* trueQualToFakeQual;
	double kappa, lambda;

	//position params
	double revIntercept = 1.0;
	double intercept = 1.0;
	double m = 0.0;

	//optimization algorithm params
	double* w;
	bool weightsInitialized;
//	double kappa_cur = -1.0, lambda_cur = -1.0;
	std::vector< std::vector<double> > QBetaQBetaP;
	bool betaQBetaPInitialized = false;

	//quality functions
	void parseBQSRQualInput(TParameters & params);
	double returnFakeError(int & trueQual);
	void setFakeQualityDistribution();
	void initializeTrueQualToFakeQualTable();
	int sampleFakeQuality();

	//position functions
	void calculateSlopeIntercept();
	double returnBetaPp(int & pos);

	//optimization algorithm functions
	void fillWeights(double & kappa_cur, double & lambda_cur);
	void fillQBetaQBetaP();
	double returnCurMean();
	double returnCurSD(double & kappa);
	double returnDelta(double & kappa);
	void simulate(short* posAddress, readLengthContainer & rl, TGenotypeMap & genoMap);

public:
	TSimulatorReadBQSR(TSimulatorReadLength* ReadLengthDist, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase);
	virtual ~TSimulatorReadBQSR(){
		if(trueQualToFakeQualTableInitialized)
			delete[] trueQualToFakeQual;
		if(weightsInitialized)
			delete[] w;
	};
};


#endif /* TSIMULATORREAD_H_ */
