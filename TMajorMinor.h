/*
 * TMajorMinor.h
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#ifndef TMAJORMINOR_H_
#define TMAJORMINOR_H_

#include <math.h>
#include "TGLF.h"
#include "TGenotypeMap.h"

class TMajorMinor{
private:
	TLog* logfile;
	TQualityMap qualiMap;

	void fillInitialEstimateOfGenotypeFrequencies(double* genoFreq, const int & numSamples, double* genotypeLikelihoods);
	void estimateGenotypeFrequencies(double* genotypeFrequencies, const int & numSamples, double* genotypeLikelihoods, const double & epsilonF);
	void calculateLogLikelihood(const int & numSamples, double* genotypeLikelihoods, double* genotypeFrequencies);
	void fillLikelihoods(const int & numSamples, uint8_t** phred, double* genotypeLikelihoods, int aa, int ab, int bb);
	void getMLEOfMajorMinor(const int & numSamples, uint8_t** phred, double* genotypeLikelihoods, const double & epsilonF);

public:
	TMajorMinor(TLog* Logfile);

	void estimateMajorMinor(TParameters & params);

};



#endif /* TMAJORMINOR_H_ */
