/*
 * TAlleleCountEstimator.h
 *
 *  Created on: Dec 9, 2018
 *      Author: phaentu
 */

#ifndef TALLELECOUNTESTIMATOR_H_
#define TALLELECOUNTESTIMATOR_H_

#include "mathFunctions.h"
#include "TPopulationLikelihoods.h"

//-------------------------------------------------
// TSiteAlleleFrequencyLikelihoods
//-------------------------------------------------
class TSiteAlleleFrequencyLikelihoods{
private:
	TQualityMap qualMap;
	int numInd_k;
	double logOf2;
	int numAlleleCounts; //from 0 to 2k
	double* log_choose_2k_j;
	double* log_alleleFrequencyLikelihoods_h;

	double protectedSumInLog(double a, double b);
	double protectedSumInLog(double a, double b, double c);
	void normalize();

	void fillLog(uint8_t* phred);
	void fillNatural(uint8_t* phred);

public:
	TSiteAlleleFrequencyLikelihoods(int numIndividuals);
	~TSiteAlleleFrequencyLikelihoods();

	void fill(uint8_t* phred);
	void print();
};

//-------------------------------------------------
// TAlleleCountEstimator
//-------------------------------------------------
class TAlleleCountEstimator{
private:
	TLog* logfile;
	TRandomGenerator* randomGenerator;

public:
	TAlleleCountEstimator(TParameters & params, TLog* Logfile);
	~TAlleleCountEstimator();

	void estimateAlleleCounts(TParameters & params);
};


#endif /* TALLELECOUNTESTIMATOR_H_ */
