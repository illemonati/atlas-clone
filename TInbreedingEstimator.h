/*
 * TInbreedingEstimator.h
 *
 *  Created on: Dec 6, 2018
 *      Author: linkv
 */

#ifndef TINBREEDINGESTIMATOR_H_
#define TINBREEDINGESTIMATOR_H_

#include "TParameters.h"
#include "TLog.h"
#include "TPopulationLikelihoods.h"
#include <limits>

//---------------------------
// alphaOrBeta
//---------------------------

class TAlphaOrBeta{
private:
	double _alphaOrBeta;
	double _logAlphaOrBeta;
public:
	std::string variableName;
	TAlphaOrBeta(std::string VariableName);
	bool initialize(std::vector<double> & p, TLog* logfile);
	void update(double & newLogValue, double & newNaturalScaleValue);
	double getLogValue();
	double getNaturalScaleValue();
};

//---------------------------
// TInbreedingEstimator
//---------------------------

class TInbreedingEstimator{
private:
	TRandomGenerator randomGenerator;

	//log
	TLog* logfile;
	std::string outname;

	//algorithm params
	int numIterations;
	double pi;
	double widthProposalKernelLogAlphaOrBeta;

	//data
	TPopulationLikelihoods likelihoods;
	unsigned long numLoci;

	//params
	double F;
	std::vector<double> p;
	TAlphaOrBeta alpha = TAlphaOrBeta("");
	TAlphaOrBeta beta = TAlphaOrBeta("");

//	void initializeAlphaBeta();
	void initParams(TRandomGenerator & randomGenerator);
	void printTrajectory(gz::ogzstream & tracefile);
	void updateF();
	void updateP(long l, TAlphaOrBeta & alpha, TAlphaOrBeta & beta);
	bool updateAlphaOrBeta(TAlphaOrBeta & alphaOrBetaToUpdate, TAlphaOrBeta & alphaOrBetaOther);
	double PGenoGivenFAndP(int & genotype, double & F, double & p);
	void oneMCMCIteration();

public:
	TInbreedingEstimator(TParameters & Parameters, TLog* Logfile);
	void runEstimation();
};



#endif /* TINBREEDINGESTIMATOR_H_ */
