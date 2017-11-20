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
	TThetaEstimatorTemporaryFile();
	TThetaEstimatorTemporaryFile(std::string filename, int numGenotypes);
	~TThetaEstimatorTemporaryFile(){
		clean();
	};

	void init(std::string filename, int numGenotypes);

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
	int maxKforPoissonPlusOne;
	double* poissonProb;
	uint8_t* numBootstrapRepsPerEntry;
	bool numBootstrapRepsPerEntryInitialized;

	bool readState;
	long curSite;
	uint8_t curRep;
	double* pointerToData;

	//tmp variables
	int g;
	double sum;
	double* P_g_oneSite;

	virtual void saveSite(TSite & site){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual void emptyStorage(){};
	void fillPoissonForBootstrap(const double lambda);
	virtual void _begin(){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual void readNext(){ throw "Not available in TThetaEstimatorData base class!"; };

public:
	TThetaEstimatorData(int NumGenotypes);
	virtual ~TThetaEstimatorData(){
		clear();
		delete[] poissonProb;
		delete[] P_g_oneSite;
		if(numBootstrapRepsPerEntryInitialized)
			delete[] numBootstrapRepsPerEntry;
	};

	void add(TSite & site);
	void clear();
	void bootstrap(TRandomGenerator & randomGenerator);
	void clearBootstrap();

	virtual bool begin();
	virtual bool next();
	virtual bool isEnd(){ throw "Not available in TThetaEstimatorData base class!"; };
	virtual double* curGenotypeLikelihoods(){ throw "Not available in TThetaEstimatorData base class!"; };

	long size(){return totNumSitesAdded;};
	long sizeWithData(){
		if(isBootstrapped)
			return numBootstrappedSites;
		return numSitesWithData;
	};

	void writeHeader(std::ofstream & out);
	void writeSize(std::ofstream & out);
	void fillBaseFreq(double* baseFreq);
	virtual void fillP_G(double* & P_G, double* & pGenotype);
	virtual double calcLogLikelihood(double* & pGenotype);
};


class TThetaEstimatorDataVector:public TThetaEstimatorData{
private:
	std::vector<double*> sites;
	std::vector<double*>::iterator siteIt;

	void saveSite(TSite & site);
	void emptyStorage();
	void readNext();

public:
	TThetaEstimatorDataVector(int NumGenotypes);
	virtual ~TThetaEstimatorDataVector(){
		clear();
	};

	void _begin();
	bool isEnd();
	double* curGenotypeLikelihoods();
	void fillP_G(double* & P_G, double* & pGenotype);
	double calcLogLikelihood(double* & pGenotype);
};


class TThetaEstimatorDataFile:public TThetaEstimatorData{
protected:
	std::string dataFileName;

	TThetaEstimatorTemporaryFile sites;

	void saveSite(TSite & site);
	void emptyStorage();
	void readNext();

public:
	TThetaEstimatorDataFile(int NumGenotypes, std::string TmpFileName);
	~TThetaEstimatorDataFile(){
		delete[] pointerToData;
		clear();
	};

	void _begin();
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

	void clear();
	void add(TSite & site);
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
