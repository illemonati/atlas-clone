/*
 * TAlleleCountEstimator.h
 *
 *  Created on: Dec 9, 2018
 *      Author: phaentu
 */

#ifndef TALLELECOUNTESTIMATOR_H_
#define TALLELECOUNTESTIMATOR_H_

#include "TPopulationLikelihoods.h"


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
