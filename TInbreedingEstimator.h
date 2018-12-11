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

class TInbreedingEstimator{
private:
	TRandomGenerator randomGenerator;

	//log
	TLog* logfile;
	std::string outname;

	//algorithm params
	int numIterations;
	double pi;

	//data
	TPopulationLikelihoods likelihoods;

	//params
	double F;
	std::vector<double> p;
	double alpha, logAlpha;
	double beta, logBeta;

	void initializeAlphaBeta();
	void initParams(TRandomGenerator & randomGenerator);
	void printTrajectory(gz::ogzstream & tracefile);
	void updateF();
	void updateP(long l);
	void updateAlphaOrBeta();

public:
	TInbreedingEstimator(TParameters & Parameters, TLog* Logfile);
	void runEstimation();
};



#endif /* TINBREEDINGESTIMATOR_H_ */
