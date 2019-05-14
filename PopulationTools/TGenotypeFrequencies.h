/*
 * TGenotypeFreqencies.h
 *
 *  Created on: Apr 30, 2019
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TGENOTYPEFREQUENCIES_H_
#define POPULATIONTOOLS_TGENOTYPEFREQUENCIES_H_

#include "TPopulationLikelihoodStorage.h"
#include "../TQualityMap.h"
#include <iostream>
#include "../TFile.h"

//------------------------------------------------
// TGenotypeFreqencies
//------------------------------------------------
class TGenotypeFrequencies{
private:

	void setFrequencenciesToZero();
	void normalize();
	void ensureAllFrequenciesAreNonZero();

public:
	double diploidFrequencies[3];
	double haploidFrequencies[2];
	double alleleFrequency;
	double MAF;
	int numDiploidSamples;
	int numHaploidSamples;

	TGenotypeFrequencies();

	void clear();
	void set(const TGenotypeFrequencies & other);
	void flip(); //flip major / minor
	bool isMonomorphic();
	void guess(TPopulationLikehoodStorage & samples);
	void guess(TSampleLikelihoods* samples, int numSamples);
	void estimate(TPopulationLikehoodStorage & samples, TQualityMap & phredToGTLMap, double epsilonF);
	void estimate(TSampleLikelihoods* samples, int numSamples, TQualityMap & phredToGTLMap, double epsilonF);
	double calculateLog10Likelihood(TPopulationLikehoodStorage & samples, TQualityMap & phredToGTLMap);
	double calculateLog10Likelihood(TSampleLikelihoods* samples, int numSamples, TQualityMap & phredToGTLMap);
	void writeDiploidFrequencies(TOutputFile & out);
	void writeHaploidFrequencies(TOutputFile & out);
};




#endif /* POPULATIONTOOLS_TGENOTYPEFREQUENCIES_H_ */
