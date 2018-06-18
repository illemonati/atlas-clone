/*
 * TThetaEstimator.h
 *
 *  Created on: Sep 25, 2017
 *      Author: phaentu
 */

#ifndef TTHETAESTIMATOR_H_
#define TTHETAESTIMATOR_H_

#include "TThetaEstimatorData.h"
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
//TThetaEstimator
//---------------------------------------------------------------
class TThetaEstimator_base{
protected:
	TLog* logfile;

	//data
	TThetaEstimatorData* data;
	bool dataInitialized;
	bool useTmpFile;
	std::string tmpFileName;
	int numGenotypes;
	TGenotypeMap genoMap;

	//initial theta
	double initalTheta;
	double initThetaSearchFactor;
	int initThetaNumSearchIterations;

	//estimation
	int minSitesWithData;
	double* baseFreq;
	double* pGenotype; //P(g|pi, theta)
	double* P_G; // see paper
	double* P_g_oneSite;
	Theta theta;
	bool extraVerbose;

	void initTmpStorage();
	void readParametersRegardingInitialSearch(TParameters & params);
	void fillPGenotype(double* & pGeno, double & expTheta);
	void fillPGenotype(double & expTheta);
	void fillP_G();
	double calcLogLikelihood();

	void findGoodStartingTheta(TThetaEstimatorData* thisData);

public:
	TThetaEstimator_base(TLog* Logfile);
	TThetaEstimator_base(TParameters & params, TLog* Logfile);

	~TThetaEstimator_base(){
		if(dataInitialized)
			delete data;
		delete[] P_G;
		delete[] pGenotype;
		delete[] baseFreq;
		delete[] P_g_oneSite;
	};

	TThetaEstimatorData* pointerToDataContainer(){ return data; };

	void fillPGenotype(double* & pGeno);

};

//---------------------------------------------------------------
//TThetaEstimator
//---------------------------------------------------------------
class TThetaEstimator:public TThetaEstimator_base{
private:
	//EM parameters
	int numIterations;
	int numThetaOnlyUpdates;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;

	void init();
	double calcFisherInfo(double* deriv_pGenotype);
	void runEMForTheta();
	void estimateConfidenceInterval();

public:
	TThetaEstimator(TParameters & params, TLog* Logfile);
	TThetaEstimator(TLog* Logfile);

	~TThetaEstimator(){};

	void clear();
	void add(TSite & site);
	long sizeWithData(){ return data->sizeWithData();};
	bool estimateTheta();
	void setTheta(double Theta);
	void setBaseFreq(TBaseFrequencies & BaseFreq);
	void writeHeader(std::ofstream & out);
	void writeThetas(std::ofstream & out);
	void writeResultsToFile(std::ofstream & out);
	void calcLikelihoodSurface(std::ofstream & out, int & steps);
	void bootstrapTheta(TRandomGenerator & randomGenerator, std::ofstream & out);
};


//---------------------------------------------------------------
//TThetaEstimatorRatio
//---------------------------------------------------------------
class TThetaEstimatorRatio:public TThetaEstimator_base{
private:
	//second data
	TThetaEstimatorData* data2;
	bool data2Initialized;

	//MCMC parameters
	int numIterations;
	int burnin;



public:
	TThetaEstimatorRatio(TParameters & params, TLog* Logfile);
	~TThetaEstimatorRatio(){
		if(data2Initialized)
			delete data2;
	};

	TThetaEstimatorData* pointerToDataContainer2(){ return data2; };

	void estimateRatio();
};

#endif /* TTHETAESTIMATOR_H_ */
