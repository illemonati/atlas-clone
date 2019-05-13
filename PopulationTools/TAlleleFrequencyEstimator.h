/*
 * TAlleleFrequencyEstimator.h
 *
 *  Created on: May 13, 2019
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_
#define POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_

/*
#include "TParameters.h"
#include "TPopulationLikelihoods.cpp"


//------------------------------------------------
//TAlleleHardyWeinbergFreqEstimator
//------------------------------------------------
class TAlleleHardyWeinbergFreqEstimator{
private:
	double f;
	int maxIter;

public:
	TAlleleHardyWeinbergFreqEstimator();
	double estimate(TPopulationLikehoodStorage storage, double epsilonF);
};

//------------------------------------------------
//TAlleleFreqEstimator
//------------------------------------------------
class TAlleleFreqEstimator{
	// about vcf-file
	std::string vcfFilename;
	bool vcfRead;
	TLog* logfile;

	// data on individuals
	TPopulationSamples samples;

	//data on loci
	long _numLoci;

public:
	TAlleleFreqEstimator(TParameters & Parameters, TLog* logfile);
    void estimateAlleleFreq(TParameters & Parameters);
};
*/


#endif /* POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_ */
