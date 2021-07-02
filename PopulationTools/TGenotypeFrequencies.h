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

using coretools::Probability;

//------------------------------------------------
// TGenotypeFreqencies
//------------------------------------------------
class TGenotypeFrequencies{
private:

	//void normalize();
	void ensureAllFrequenciesAreNonZero();

public:
	std::array<Probability, 3> diploidFrequencies;
	std::array<Probability, 2> haploidFrequencies;
	Probability alleleFrequency;
	Probability MAF;
	uint32_t numDiploidSamples;
	uint32_t numHaploidSamples;

	TGenotypeFrequencies();

	void clear();
	void set(const TGenotypeFrequencies & other);
	void flip(); //flip major / minor
	bool isMonomorphic();
	void guess(TSampleLikelihoods* samples, int numSamples);
	void estimate(TPopulationLikehoodLocus & samples, const double & epsilonF);
	void estimate(TSampleLikelihoods* samples, const uint32_t & numSamples, const double & epsilonF);
	coretools::Log10Probability calculateLog10Likelihood(TPopulationLikehoodLocus & samples);
	coretools::Log10Probability calculateLog10Likelihood(TSampleLikelihoods* samples, const uint32_t & numSamples);
	void writeDiploidFrequencies(coretools::TOutputFile & out);
	void writeHaploidFrequencies(coretools::TOutputFile & out);
	uint32_t numHaploid(){ return numHaploidSamples; };
	uint32_t numDiploid(){ return numDiploidSamples; };
};

}; //end namespace

#endif /* POPULATIONTOOLS_TGENOTYPEFREQUENCIES_H_ */
