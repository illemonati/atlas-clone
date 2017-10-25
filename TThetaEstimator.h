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
//TThetaEstimator
//---------------------------------------------------------------
class TThetaEstimator{
private:
	TLog* logfile;

	//data
	std::vector<double*> sites;
	long numSitesCoveredTwiceOrMore;
	long totNumSitesAdded;
	long numSitesWithData;
	double cumulativeDepth;

	//EM parameters
	int numIterations;
	int numThetaOnlyUpdates;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
	double initalTheta;
	double initThetaSearchFactor;
	int initThetaNumSearchIterations;

	//estimation
	int numGenotypes;
	double* P_G; // see paper
	double* pGenotype; //P(g|pi, theta)
	double* baseFreq;
	TGenotypeMap genoMap;
	TBaseFrequencies initialBaseFreq;
	Theta theta;

	//tmp variables
	int g;
	std::vector<double*>::iterator siteIt;
	std::vector<double*>::reverse_iterator revSiteIt;
	double* doublePointer;
	double sum;
	double P_g_oneSite[10];

	void init();
	void fillPGenotype(double* & pGeno, double & expTheta);
	void fillPGenotype(double & expTheta);
	void fillP_G(std::vector<double*> & theseSites);
	double calcLogLikelihood(std::vector<double*> & theseSites);
	void findGoodStartingTheta(std::vector<double*> & theseSites);
	void runEMForTheta(std::vector<double*> & theseSites);
	void estimateConfidenceInterval(std::vector<double*> & theseSites);
	bool estimateTheta(std::vector<double*> & theseSites);

public:
	TThetaEstimator(TParameters & params, TLog* Logfile);
	TThetaEstimator(TLog* Logfile);

	~TThetaEstimator(){
		clear();
		delete[] P_G;
		delete[] pGenotype;
		delete[] baseFreq;
	};

	void clear(){
		for(siteIt=sites.begin(); siteIt != sites.end(); ++siteIt)
			delete[] (*siteIt);
		sites.clear();
		initialBaseFreq.clear();
		numSitesCoveredTwiceOrMore = 0;
		totNumSitesAdded = 0;
		numSitesWithData = 0;
		cumulativeDepth = 0.0;
	};
	void add(TSite & site);
	long size(){ return sites.size(); };
	void fillPGenotype(double* & pGeno);
	bool estimateTheta();
	void setTheta(double Theta);
	void setBaseFreq(TBaseFrequencies & BaseFreq);
	void writeHeader(std::ofstream & out);
	void writeResultsToFile(std::ofstream & out);
	void calcLikelihoodSurface(std::ofstream & out, int & steps);
	void bootstrapTheta(TRandomGenerator & randomGenerator, std::ofstream & out);
};


#endif /* TTHETAESTIMATOR_H_ */
