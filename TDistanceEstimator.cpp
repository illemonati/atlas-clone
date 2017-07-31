/*
 * TDistanceCalculator.cpp
 *
 *  Created on: Jul 25, 2017
 *      Author: phaentu
 */

#include "TDistanceEstimator.h"

void TDistanceEstimator::calcDistance(TParameters & params){

	//test first to parse GLF files
	std::string glf = params.getParameterString("glf");
	TGlfReader reader(glf);

	//print file
	reader.printToEnd();
}


