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
#include "TGenotypeMap.h"

class TDistanceEstimate{
private:
	TGenotypeMap genoMap;
	TPhredToLikelihood phredToLik;

public:
	TBaseFrequencies pi, old_pi;
	double* phi;
	double* old_phi;
	double LL;
	double old_LL;

	TDistanceEstimate();
	void guessPi(int** genoQual1, int** genoQual2, long numSites);

};

class TDistanceEstimator{
private:
	TLog* logfile;
	TPhredToLikelihood phredToLik;

	void initializePhredToLikelihoodTable();
	void deletePhredToLikelihoodTable();
	void estimateDistanceInWindows(std::string filename, TGlfReader & g1, TGlfReader & g2, long windowLen);

	void writeDistanceEstimates(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd, int & numsitesWithData);
	void writeDistanceEstimatesNoData(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd);

public:
	TDistanceEstimator(TLog* Logfile);
	~TDistanceEstimator(){};

	void printGLF(TParameters & params);
	void estimateDistances(TParameters & params);

};


#endif /* TDISTANCEESTIMATOR_H_ */
