/*
 * TGenotypeFreqencies.cpp
 *
 *  Created on: Apr 30, 2019
 *      Author: wegmannd
 */

#include "TGenotypeFrequencies.h"
#include "debugtools.h"

namespace PopulationTools{

////////////////////////////////////////////////////////////////////////////////////////////////
// TGenotypeFreqencies                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////
TGenotypeFrequencies::TGenotypeFrequencies(){
	clear();
};

void TGenotypeFrequencies::clear(){
	numDiploidSamples = 0;
	numHaploidSamples = 0;
};

void TGenotypeFrequencies::set(const TGenotypeFrequencies & other){
	diploidFrequencies[0] = other.diploidFrequencies[0];
	diploidFrequencies[1] = other.diploidFrequencies[1];
	diploidFrequencies[2] = other.diploidFrequencies[2];

	haploidFrequencies[0] = other.haploidFrequencies[0];
	haploidFrequencies[1] = other.haploidFrequencies[1];

	alleleFrequency = other.alleleFrequency;
	MAF = other.MAF;

	numDiploidSamples = other.numDiploidSamples;
	numHaploidSamples = other.numHaploidSamples;
};

void TGenotypeFrequencies::flip(){
	//flip major / minor
	double tmp = haploidFrequencies[0];
	haploidFrequencies[0] = haploidFrequencies[1];
	haploidFrequencies[1] = tmp;

	tmp = diploidFrequencies[0];
	diploidFrequencies[1] = diploidFrequencies[0];
	diploidFrequencies[0] = diploidFrequencies[1];

	alleleFrequency = 1.0 - alleleFrequency;
};

bool TGenotypeFrequencies::isMonomorphic(){
	bool isMonoDiploid = numDiploidSamples > 0 && (diploidFrequencies[0] == 1.0 || diploidFrequencies[1] == 0.0);
	bool isMonoHaploid = numHaploidSamples > 0 && (haploidFrequencies[0] == 1.0 || haploidFrequencies[1] == 0.0);

	return isMonoDiploid && isMonoHaploid;
};

void TGenotypeFrequencies::ensureAllFrequenciesAreNonZero(){
	//diploid
	if(numDiploidSamples > 0){
		if(diploidFrequencies[0] == 0.0)
			diploidFrequencies[0] = 0.000001;
		if(diploidFrequencies[1] == 0.0)
			diploidFrequencies[1] = 0.000001;
		if(diploidFrequencies[2] == 0.0)
			diploidFrequencies[2] = 0.000001;

		//renormalize
		double sum = (double) diploidFrequencies[0] + (double) diploidFrequencies[1] + (double) diploidFrequencies[2];
		diploidFrequencies[0] /= sum;
		diploidFrequencies[1] /= sum;
		diploidFrequencies[2] /= sum;
	}

	//haploid
	if(numHaploidSamples > 0){
		if(haploidFrequencies[0] == 0.0)
			haploidFrequencies[0] = 0.000001;
		if(haploidFrequencies[1] == 0.0)
			haploidFrequencies[1] = 0.000001;

		//renormalize
		double sum = (double) haploidFrequencies[0] + (double) haploidFrequencies[1];
		haploidFrequencies[0] /= sum;
		haploidFrequencies[1] /= sum;
	}
};

void TGenotypeFrequencies::estimate(TPopulationLikehoodLocus & samples, const double  &epsilonF){
	estimate(samples.samples(), samples.numSamples(), epsilonF);
};

coretools::Log10Probability TGenotypeFrequencies::calculateLog10Likelihood(TPopulationLikehoodLocus & samples){
	return calculateLog10Likelihood(samples.samples(), samples.numSamples());
};

void TGenotypeFrequencies::writeDiploidFrequencies(coretools::TOutputFile & out){
	if(numDiploidSamples > 0){
		out << diploidFrequencies[0] << diploidFrequencies[1] << diploidFrequencies[2];
	} else {
		out << "-" << "-" << "-";
	}
};

void TGenotypeFrequencies::writeHaploidFrequencies(coretools::TOutputFile & out){
	if(numHaploidSamples > 0){
		out << haploidFrequencies[0] << haploidFrequencies[1];
	} else {
		out << "-" << "-";
	}
};

}; //end namespace
