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
	bool dePhredTableInitialized = false;
	double* dePhredTable;
	double meanQual, sdQual;
	int maxQual;
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

class TSimulatorReadBQSRPos:public TSimulatorRead{
private:

protected:
	TSimulatorReadLength* readLengthDist;
	int maxPos;

	//quality params
	int phi1;
	double phi2;
	bool fakeQualToTrueQualTableInitialized = false;
	double* fakeQualToTrueQual;
	double kappa, lambda;

	//position params
	double revIntercept;
	double intercept;
	double m;

	//optimization algorithm params
	int minQual = 0;
	double* w;
	bool weightsInitialized;
	double kappa_cur = -1.0, lambda_cur = -1.0;
	double delta, delta_old;
	std::vector< std::vector<double> > QBetaQBetaP;
	bool BetaQBetaPInitialized = false;

	//quality functions
	void parseBQSRQualInput(TParameters & params);
	int returnTrueQual(int & fakeQual);
	void setFakeQualityDistribution();
	void initializeFakeQualToTrueQualTable();
	int sampleFakeQuality();

	//position functions
	void calculateSlopeIntercept();
	double returnBetaPp(int pos);

	//optimization algorithm functions
	void fillWeights(double & kappa_cur, double & lambda_cur);
	void fillQBetaQBetaP();
	double returnCurKappa();
	double returnCurLambda();
	double returnDelta();
	void simulate(short* posAddress, readLengthContainer & rl, TGenotypeMap & genoMap);

public:
	TSimulatorReadBQSRPos(TSimulatorReadLength* ReadLengthDist, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase);
	virtual ~TSimulatorReadBQSRPos(){
		if(fakeQualToTrueQualTableInitialized)
			delete[] fakeQualToTrueQual;
		if(weightsInitialized)
			delete[] w;
	};
};

//-------------------------------
// BQSR transformation qual
//-------------------------------

class TSimulatorReadBQSRQual:public TSimulatorReadBQSRPos{
public:
	TSimulatorReadBQSRQual(TSimulatorReadLength* ReadLengthDist, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase);
	virtual ~TSimulatorReadBQSRQual(){};
};


#endif /* TSIMULATORREAD_H_ */
