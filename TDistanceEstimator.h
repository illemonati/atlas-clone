/*
 * TDistanceCalculator.h
 *
 *  Created on: Jul 25, 2017
 *      Author: phaentu
 */

#ifndef TDISTANCEESTIMATOR_H_
#define TDISTANCEESTIMATOR_H_

#include "TParameters.h"
#include "TGLF.h"
#include <math.h>

class TDistanceEstimator{
private:
	TLog* logfile;
	double* phredToLikelihood;
	bool phredToLikelihoodInitialized;

	void initializePhredToLikelihoodTable();
	void deletePhredToLikelihoodTable();
	void estimateDistanceInWindows(std::string filename, TGlfReader & g1, TGlfReader & g2, long windowLen);

	void writeDistanceEstimates(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd, int & numsitesWithData);
	void writeDistanceEstimatesNoData(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd);

public:
	TDistanceEstimator(TLog* Logfile);
	~TDistanceEstimator(){
		deletePhredToLikelihoodTable();
	};

	void printGLF(TParameters & params);
	void estimateDistances(TParameters & params);

};


#endif /* TDISTANCEESTIMATOR_H_ */
