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
//TThetaEstimatorDataStorage
//---------------------------------------------------------------
class TThetaEstimatorDataStorage{
private:
	std::vector<double*> sites;
	std::vector<double*>* useTheseSites;
	std::vector<double*>::iterator siteIt;
	std::vector<double*>::reverse_iterator revSiteIt;
	double* doublePointer;

	std::vector<double*> bootstrappedSites;

protected:
	int numGenotypes;
	long numSitesCoveredTwiceOrMore;
	long totNumSitesAdded;
	long numSitesWithData;
	double cumulativeDepth;

	TBaseFrequencies initialBaseFreq;
	bool isBootstrapped;
	long numBootstrappedSites;

	//tmp variables
	int g;
	double sum;
	double* P_g_oneSite;

public:
	TThetaEstimatorDataStorage(int NumGenotypes);
	virtual ~TThetaEstimatorDataStorage(){
		clear();
		delete[] P_g_oneSite;
	}

	virtual void add(TSite & site);
	virtual void bootstrap(TRandomGenerator & randomGenerator);
	virtual void clearBootstrap();
	virtual void clear();

	long size(){return totNumSitesAdded;};
	long sizeWithData(){return numSitesWithData;};

	void writeHeader(std::ofstream & out);
	void writeSize(std::ofstream & out);

	void fillBaseFreq(double* baseFreq);
	void fillP_G(double* P_G, double* pGenotype);
	double calcLogLikelihood(double* pGenotype);
	double calcFisherInfo(double* pGenotype, double* deriv_pGenotype);
};

/*
class TThetaEstimatorDataStorageFile:public TThetaEstimatorDataStorage{
private:
	bool dataStoredInFile;
	std::string dataFileName;
	gzFile fp;
	bool dataFileOpen;

public:
	TThetaEstimatorDataStorageFile();
	~TThetaEstimatorDataStorageFile();
	void add(TSite & site);
	void clear();
};
*/

//---------------------------------------------------------------
//TThetaEstimator
//---------------------------------------------------------------
class TThetaEstimator{
private:
	TLog* logfile;

	TThetaEstimatorDataStorage* data;

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
	double P_g_oneSite[10];

	void init();
	void fillPGenotype(double* & pGeno, double & expTheta);
	void fillPGenotype(double & expTheta);
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
