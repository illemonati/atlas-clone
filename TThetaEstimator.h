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
	double theta, expTheta, logTheta, thetaConfidence, LL;
	double* baseFreq;

	Theta(){
		theta = 0.0;
		thetaConfidence = 0.0;
		expTheta = 0.0;
		logTheta = -9e100;
		LL = -9e100;
		baseFreq = new double[4];
	};

	~Theta(){
		delete[] baseFreq;
	};

	void setTheta(double val){
		theta = val;
		expTheta = exp(-theta);
		logTheta = log(theta);
	};

	void setExpTheta(double val){
		expTheta = val;
		theta = -log(val);
		logTheta = log(theta);
	};

	void setLogTheta(double val){
		logTheta = val;
		theta = exp(val);
		expTheta = exp(-theta);
	};

	void setLogTheta(double & val, double & newLL){
		logTheta = val;
		theta = exp(val);
		expTheta = exp(-theta);
		LL = newLL;
	};

	std::string getBaseFrequencyString(){
		return "Pi(A) = " + toString(baseFreq[0]) + ", Pi(C) = " + toString(baseFreq[1]) + ", Pi(G) = " + toString(baseFreq[2]) + ", Pi(T) = " + toString(baseFreq[3]);
	};
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
	double* pGenotype; //P(g|pi, theta)
	Theta theta;
	bool extraVerbose;

	void initTmpStorage();
	void readParametersRegardingInitialSearch(TParameters & params);
	void fillPGenotype(double* & pGeno, const double & expTheta, const double* baseFrequencies);
	void fillPGenotype(double* & pGeno, const Theta & thisTheta);

	void findGoodStartingTheta(TThetaEstimatorData* thisData, Theta & thisTheta, std::string tag);

public:
	TThetaEstimator_base(TLog* Logfile);
	TThetaEstimator_base(TParameters & params, TLog* Logfile);

	virtual ~TThetaEstimator_base(){
		if(dataInitialized)
			delete data;
		delete[] pGenotype;
	};

	TThetaEstimatorData* pointerToDataContainer(){ return data; };

	void fillPGenotype(double* & pGeno){ fillPGenotype(pGeno, theta); };

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

	//tmp vectors
	double* P_G; // see paper
	double* P_g_oneSite;

	void initAdditionalTmpStorage();
	void fillP_G();
	double calcFisherInfo(double* _pGenotype, double* deriv_pGenotype);
	void runEMForTheta();
	void estimateConfidenceInterval();

public:
	TThetaEstimator(TParameters & params, TLog* Logfile);
	TThetaEstimator(TLog* Logfile);

	virtual ~TThetaEstimator(){
		delete[] P_G;
		delete[] P_g_oneSite;
	};

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
	Theta theta2;
	double* tmpBaseFreq;

	//MCMC parameters
	double phiPriorMean;
	double phiPriorVar;
	double phiPriorOneOverTwoVar;
	double sdProposalKernelTheta1;
	double sdProposalKernelTheta2;
	double sdProposalKernelBaseFreq1;
	double sdProposalKernelBaseFreq2;
	int numIterations;
	int burnin;
	int thinning;
	int numAcceptedTheta1;
	int numAcceptedTheta2;
	int numAcceptedBaseFreq1;
	int numAcceptedBaseFreq2;

	void initAdditionalTmpStorage();
	void clearCounters();
	void concludeAcceptanceRate(const int & numAccepted, const int & length, std::string name);
	void concludeAcceptanceRateUpdateProposal(const int & numAccepted, const int & length, double & sd, std::string name);
	void concludeAcceptanceRates(const int & length);
	void concludeAcceptanceRatesUpdateProposal(const int & length);
	bool updateTheta(TThetaEstimatorData* thisData, Theta & thisTheta, double otherLogThetaMean, const double & thisSdProposalKernel, TRandomGenerator & randomGenerator);
	bool updateBaseFrequencies(TThetaEstimatorData* thisData, Theta & thisTheta, const double & thisSdProposalKernel, TRandomGenerator & randomGenerator);
	void oneMCMCIteration(TRandomGenerator & randomGenerator);

public:
	TThetaEstimatorRatio(TParameters & params, TLog* Logfile);
	~TThetaEstimatorRatio(){
		if(data2Initialized)
			delete data2;
		delete[] tmpBaseFreq;
	};

	TThetaEstimatorData* pointerToDataContainer2(){ return data2; };

	void estimateRatio(TRandomGenerator & randomGenerator, std::string ouputName);


};

#endif /* TTHETAESTIMATOR_H_ */
