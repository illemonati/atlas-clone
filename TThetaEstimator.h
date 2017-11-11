/*
 * TThetaEstimator.h
 *
 *  Created on: Sep 25, 2017
 *      Author: phaentu
 */

#ifndef TTHETAESTIMATOR_H_
#define TTHETAESTIMATOR_H_

#include "TSite.h"
#include "TRandomGenerator.h"
#include <stdio.h>

//---------------------------------------------------------------
//Theta
//---------------------------------------------------------------
struct Theta{
	double theta, expTheta, thetaConfidence, LL;

	Theta(){
		theta = 0.0;
		thetaConfidence = 0.0;
		expTheta = 0.0;
		LL = -9e100;
	};

	void setTheta(double val){
		theta = val;
		expTheta = exp(-theta);
	}
};

//---------------------------------------------------------------
//TThetaEstimatorTemporaryFile
//---------------------------------------------------------------
class TThetaEstimatorTemporaryFile{
private:
	std::string filename;
	gzFile fp;
	int sizeOfData;

	bool isOpen;
	bool isOpenForWriting;
	bool isOpenForReading;
	bool wasWritten;

public:
	TThetaEstimatorTemporaryFile(std::string filename, int numGenotypes);
	~TThetaEstimatorTemporaryFile(){
		clean();
	};

	void openForWriting();
	void openForReading();
	void close();
	void clean();
	bool isEOF();

	void save(double* data);
	bool read(double* data);
};


//---------------------------------------------------------------
//TThetaEstimatorData
//---------------------------------------------------------------
class TThetaEstimatorData{
protected:
	//counters
	int numGenotypes;
	long numSitesCoveredTwiceOrMore;
	long totNumSitesAdded;
	long numSitesWithData;
	double cumulativeDepth;

	TBaseFrequencies initialBaseFreq;
	bool isBootstrapped;
	long numBootstrappedSites;

	double* pointerToData;

	virtual void saveSite(TSite & site){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual void emptyStorage(){};

public:
	TThetaEstimatorData(int NumGenotypes);
	virtual ~TThetaEstimatorData(){
		clear();
	};

	void add(TSite & site);
	void clear();
	virtual void bootstrap(TRandomGenerator & randomGenerator){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual void clearBootstrap(){};

	virtual bool begin(){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual bool next(){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual bool isEnd(){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual double* curGenotypeLikelihoods(){ throw "Not available in TThetaEstimatorData base class!"; };

	long size(){return totNumSitesAdded;};
	long sizeWithData(){return numSitesWithData;};

	void writeHeader(std::ofstream & out);
	void writeSize(std::ofstream & out);

	void fillBaseFreq(double* baseFreq);
};


class TThetaEstimatorDataVector:public TThetaEstimatorData{
private:
	std::vector<double*> sites;
	std::vector<double*> bootstrappedSites;
	std::vector<double*>* useTheseSites;
	std::vector<double*>::iterator siteIt;

	void saveSite(TSite & site);
	void emptyStorage();

public:
	TThetaEstimatorDataVector(int NumGenotypes);
	virtual ~TThetaEstimatorDataVector(){
		clear();
	};

	void bootstrap(TRandomGenerator & randomGenerator);
	void clearBootstrap();

	bool begin();
	bool next();
	bool isEnd();
	double* curGenotypeLikelihoods();
};


class TThetaEstimatorDataFile:public TThetaEstimatorData{
protected:
	std::string dataFileName;
	std::string bootstrapFileName;

	TThetaEstimatorTemporaryFile* sites;
	TThetaEstimatorTemporaryFile* bootstrappedSites;
	TThetaEstimatorTemporaryFile* useTheseSites;

	int maxKforPoissonPlusOne;
	double* poissonProb;

	void saveSite(TSite & site);
	void emptyStorage();

public:
	TThetaEstimatorDataFile(int NumGenotypes, std::string TmpFileName);
	~TThetaEstimatorDataFile(){
		delete[] pointerToData;
		delete[] poissonProb;
		clear();

		delete sites;
		delete bootstrappedSites;
	};

	void fillPoissonForBootstrap(const long toPick, const long remaining);
	void bootstrap(TRandomGenerator & randomGenerator);
	void clearBootstrap();

	bool begin();
	bool next();
	bool isEnd();
	double* curGenotypeLikelihoods();
};


//---------------------------------------------------------------
//TThetaEstimator
//---------------------------------------------------------------
class TThetaEstimator{
private:
	TLog* logfile;

	TThetaEstimatorData* data;

	//EM parameters
	int numIterations;
	int numThetaOnlyUpdates;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
	double initalTheta;
	double initThetaSearchFactor;
	int initThetaNumSearchIterations;
	bool extraVerbose;

	//estimation
	int numGenotypes;
	double* P_G; // see paper
	double* pGenotype; //P(g|pi, theta)
	double* baseFreq;
	TGenotypeMap genoMap;
	Theta theta;

	//tmp variables
	int g;
	double sum;
	double* P_g_oneSite;

	void init();
	void fillPGenotype(double* & pGeno, double & expTheta);
	void fillPGenotype(double & expTheta);
	void fillP_G();
	double calcLogLikelihood();
	double calcFisherInfo(double* deriv_pGenotype);
	void findGoodStartingTheta();
	void runEMForTheta();
	void estimateConfidenceInterval();

public:
	TThetaEstimator(TParameters & params, TLog* Logfile);
	TThetaEstimator(TLog* Logfile);

	~TThetaEstimator(){
		delete[] P_G;
		delete[] pGenotype;
		delete[] baseFreq;
		delete[] P_g_oneSite;
		delete data;
	};

	void add(TSite & site);
	//long size(){ return data->size();};
	long sizeWithData(){ return data->sizeWithData();};
	void fillPGenotype(double* & pGeno);
	bool estimateTheta();
	void setTheta(double Theta);
	void setBaseFreq(TBaseFrequencies & BaseFreq);
	void writeHeader(std::ofstream & out);
	void writeThetas(std::ofstream & out);
	void writeResultsToFile(std::ofstream & out);
	void calcLikelihoodSurface(std::ofstream & out, int & steps);
	void bootstrapTheta(TRandomGenerator & randomGenerator, std::ofstream & out);
};


#endif /* TTHETAESTIMATOR_H_ */
