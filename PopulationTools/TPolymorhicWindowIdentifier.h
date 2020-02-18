/*
 * TPolymorhhicWindowIdentifier.h
 *
 *  Created on: Feb 14, 2020
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TPOLYMORHICWINDOWIDENTIFIER_H_
#define POPULATIONTOOLS_TPOLYMORHICWINDOWIDENTIFIER_H_

#include "TPopulationLikelihoods.h"

class TPolymorhicWindowIdentifier{
private:
	TLog* logfile;
	TGlfConverter glfConverter;

	double _calcQualityPolymorphic(const TPopulationLikehoodWindow & window, const uint32_t & i);

public:
	TPolymorhicWindowIdentifier(TParameters & Parameters, TLog* logfile);
	void identifyPolymorphicWindows(TParameters & Parameters, TRandomGenerator* randomGenerator);
};



#endif /* POPULATIONTOOLS_TPOLYMORHICWINDOWIDENTIFIER_H_ */
