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

class TDistanceEstimator{
private:

public:
	TDistanceEstimator(){}

	void calcDistance(TParameters & params);

};


#endif /* TDISTANCEESTIMATOR_H_ */
