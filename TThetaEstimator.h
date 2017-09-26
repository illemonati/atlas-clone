/*
 * TThetaEstimator.h
 *
 *  Created on: Sep 25, 2017
 *      Author: phaentu
 */

#ifndef TTHETAESTIMATOR_H_
#define TTHETAESTIMATOR_H_

#include "TSite.h"


//---------------------------------------------------------------
//EMParameters
//---------------------------------------------------------------
struct EMParameters{
	int numIterations;
	int numThetaOnlyUpdates;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
	double initalTheta;
	double initThetaSearchFactor;
	int initThetaNumSearchIterations;

	EMParameters();
	EMParameters(TParameters & params, TLog* logfile);
	~EMParameters(){};

	void report(TLog* logfile);
};


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
class TThetaEstimator{
private:
	std::vector<double*> sites;
	long numSitesCoveredTwiceOrMore;
	long totNumSitesAdded;
	long numSitesWithData;
	int numGenotypes;
	double P_G[10]; // see paper
	double pGenotype[10]; //P(g|pi, theta)
	double baseFreq[4];
	TGenotypeMap genoMap;
	TBaseFrequencies initialBaseFreq;
	Theta theta;

	//tmp variables
	int g;
	std::vector<double>::iterator siteIt;
	double* doublePointer;
	double sum;
	double P_g_oneSite[10];

	void fillPGenotype(double & expTheta);
	void fillP_G();
	double calcLogLikelihood();
	void findGoodStartingTheta(EMParameters & EMParams);
	void runEMForTheta(EMParameters & EMParams, long & lengthWithData);
	void estimateConfidenceInterval();

public:
	TThetaEstimator(){
		numGenotypes = 10;
		numSitesCoveredTwiceOrMore = 0;
		totNumSitesAdded = 0;
		numSitesWithData = 0;

		//tmp stuff
		g = 0;
		doublePointer = NULL;
		sum = 0.0;
	};

	void clear(){
		for(siteIt=sites.begin(); siteIt != sites.end(); ++siteIt)
			delete[] *siteIt;
		sites.clear();
		initialBaseFreq.clear();
		numSitesCoveredTwiceOrMore = 0;
		totNumSitesAdded = 0;
		numSitesWithData = 0;
	};
	void add(TSite & site);

	void calcLikelihoodSurface(std::ofstream & out, int & steps);
};


#endif /* TTHETAESTIMATOR_H_ */
