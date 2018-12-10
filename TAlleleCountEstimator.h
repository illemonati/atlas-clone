/*
 * TAlleleCountEstimator.h
 *
 *  Created on: Dec 9, 2018
 *      Author: phaentu
 */

#ifndef TALLELECOUNTESTIMATOR_H_
#define TALLELECOUNTESTIMATOR_H_

#include "TPopulationLikelihoods.h"


//-------------------------------------------------
// TSiteAlleleFrequencyLikelihoods
//-------------------------------------------------
class TSiteAlleleFrequencyLikelihoods{
private:
	int numInd_k;
	int numAlleleCounts; //from 0 to 2k
	double* alleleFrequencyLikelihoods_h;

public:
	TSiteAlleleFrequencyLikelihoods(int numIndividuals);
	~TSiteAlleleFrequencyLikelihoods();
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
