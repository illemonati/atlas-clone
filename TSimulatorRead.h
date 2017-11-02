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
	double dePhredAscii(double x);
	double dePhred(double x);
	int phred(double x);

	//virtual functions
	virtual int returnBamQual(int qual, int pos, BaseContext baseContext, int maxQual);

	//qual params
	double meanQual, sdQual;
	int maxQual;
	bool qualToErroTableInitialized = false;
	double* qualToErroTable;

	//general functions
	void setTrueQualityDistribution(double mean, double sd);
	virtual int sampleTrueQuality();
	void initializeQualToErrorTable();
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

	//derived functions
	int returnQual(int & qual, int & pos, int & context);

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
	double revIntercept;
	double intercept;
	double m;
	TSimulatorReadLength* readLengthDist;

	double phi1, phi2;
	bool fakeQualToTrueQualTableInitialized = false;
	double* fakeQualToTrueQual;
	double lambda, kappa;
	void parseBQSRQualInput(TParameters & params);
	int returnTrueQual(int & fakeQual);
	void setFakeQualityDistribution();
	int sampleFakeQuality();
	void initializeFakeQualToTrueQualTable();
	void calculateSlopeIntercept();
	double returnBetaPp(int pos);
	void simulate(short* posAddress, readLengthContainer & rl, TGenotypeMap & genoMap);

public:
	TSimulatorReadBQSRPos(TSimulatorReadLength* ReadLengthDist, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase);
	virtual ~TSimulatorReadBQSRPos(){};
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
