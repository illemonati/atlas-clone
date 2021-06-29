/*
 * TGenotypeFreqencies.h
 *
 *  Created on: Apr 30, 2019
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TGENOTYPEFREQUENCIES_H_
#define POPULATIONTOOLS_TGENOTYPEFREQUENCIES_H_

#include <TPopulationLikelihoodLocus.h>
#include <iostream>
#include "TFile.h"
#include "../GLF/TGLF.h"

namespace PopulationTools{

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
	void guess(const TSampleLikelihoods* samples, int numSamples);
	void estimate(TPopulationLikehoodLocus & samples, double epsilonF);
	void estimate(TSampleLikelihoods* samples, int numSamples, double epsilonF);
	coretools::Log10Probability calculateLog10Likelihood(TPopulationLikehoodLocus & samples);
	coretools::Log10Probability calculateLog10Likelihood(TSampleLikelihoods* samples, const uint32_t & numSamples);
	void writeDiploidFrequencies(coretools::TOutputFile & out);
	void writeHaploidFrequencies(coretools::TOutputFile & out);
	int numHaploid(){ return numHaploidSamples; };
	int numDiploid(){ return numDiploidSamples; };
};

}; //end namespace

#endif /* POPULATIONTOOLS_TGENOTYPEFREQUENCIES_H_ */
