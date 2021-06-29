/*
 * TGenotypeFreqencies.cpp
 *
 *  Created on: Apr 30, 2019
 *      Author: wegmannd
 */

#include "TGenotypeFrequencies.h"

namespace PopulationTools{

////////////////////////////////////////////////////////////////////////////////////////////////
// TGenotypeFreqencies                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////
TGenotypeFrequencies::TGenotypeFrequencies(){
	clear();
};

void TGenotypeFrequencies::clear(){
	setFrequencenciesToZero();

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

void TGenotypeFrequencies::setFrequencenciesToZero(){
	diploidFrequencies[0] = 0.0;
	diploidFrequencies[1] = 0.0;
	diploidFrequencies[2] = 0.0;

	haploidFrequencies[0] = 0.0;
	haploidFrequencies[1] = 0.0;

	alleleFrequency = 0.0;
	MAF = 0.0;
};

void TGenotypeFrequencies::normalize(){
	if(numDiploidSamples > 0){
		double sum = diploidFrequencies[0] + diploidFrequencies[1] + diploidFrequencies[2];
		diploidFrequencies[0] /= sum;
		diploidFrequencies[1] /= sum;
		diploidFrequencies[2] /= sum;
	}

	if(numHaploidSamples > 0){
		double sum = haploidFrequencies[0] + haploidFrequencies[1];
		haploidFrequencies[0] /= sum;
		haploidFrequencies[1] /= sum;
	}
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
	}

	//haploid
	if(numHaploidSamples > 0){
		if(haploidFrequencies[0] == 0.0)
			haploidFrequencies[0] = 0.000001;
		if(haploidFrequencies[1] == 0.0)
			haploidFrequencies[1] = 0.000001;
	}

	//normalize
	normalize();
};

void TGenotypeFrequencies::guess(const TSampleLikelihoods* samples, int numSamples){
	//calculate by using MLE genotype for each individual
	clear();

	for(int i = 0; i < numSamples; i++){
		if(!samples[i].isMissing){
			if(samples[i].isHaploid){
				if(samples[i].glfLikelihood_1 < samples[i].glfLikelihood_0) haploidFrequencies[1] += 1.0;
				else haploidFrequencies[0] += 1.0;
				++numHaploidSamples;
			} else {
				if(samples[i].glfLikelihood_1 < samples[i].glfLikelihood_0){
					if(samples[i].glfLikelihood_2 < samples[i].glfLikelihood_1) diploidFrequencies[2] += 1.0;
					else diploidFrequencies[1] += 1.0;
				} else {
					if(samples[i].glfLikelihood_2 < samples[i].glfLikelihood_0) diploidFrequencies[2] += 1.0;
					else diploidFrequencies[0] += 1.0;
				}
				++numDiploidSamples;
			}
		}
	}

	//normalize
	normalize();
};

void TGenotypeFrequencies::estimate(TPopulationLikehoodLocus & samples, double epsilonF){
	estimate(samples.samples(), samples.numSamples(), epsilonF);
};

void TGenotypeFrequencies::estimate(TSampleLikelihoods* samples, int numSamples, double epsilonF){
	//prepare variables
	double weights[3];
	double diploidFrequencies_old[3];
	double haploidFrequencies_old[2];

	//estimate initial frequencies from MLEs
	guess(samples, numSamples);
	ensureAllFrequenciesAreNonZero();

	//run EM for max 1000 steps
	double maxF = 10000.0;
	int s = 0;

	while(maxF > epsilonF && s < 1000){
		++s;

		//set genofreq
		diploidFrequencies_old[0] = diploidFrequencies[0];
		diploidFrequencies_old[1] = diploidFrequencies[1];
		diploidFrequencies_old[2] = diploidFrequencies[2];
		haploidFrequencies_old[0] = haploidFrequencies[0];
		haploidFrequencies_old[1] = haploidFrequencies[1];
		setFrequencenciesToZero();

		//estimate new genotype frequencies
		for(int i = 0; i < numSamples; i++){
			if(!samples[i].isMissing){
				if(samples[i].isHaploid){
					weights[0] = (double) (coretools::Probability) samples[i].glfLikelihood_0 * haploidFrequencies_old[0];
					weights[1] = (double) (coretools::Probability)samples[i].glfLikelihood_1 * haploidFrequencies_old[1];

					double sum = weights[0] + weights[1];

					haploidFrequencies[0] += weights[0] / sum;
				} else {
					weights[0] = (double) (coretools::Probability) samples[i].glfLikelihood_0 * diploidFrequencies_old[0];
					weights[1] = (double) (coretools::Probability) samples[i].glfLikelihood_1 * diploidFrequencies_old[1];
					weights[2] = (double) (coretools::Probability) samples[i].glfLikelihood_2 * diploidFrequencies_old[2];

					double sum = weights[0] + weights[1] + weights[2];

					diploidFrequencies[0] += weights[0] / sum;
					diploidFrequencies[2] += weights[2] / sum;
				}
			}
		}

		if(numDiploidSamples > 0){
			diploidFrequencies[0] /= (double) numDiploidSamples;
			diploidFrequencies[2] /= (double) numDiploidSamples;
			diploidFrequencies[1] = 1.0 - diploidFrequencies[0] - diploidFrequencies[2];
		}

		if(numHaploidSamples > 0){
			haploidFrequencies[0] /= (double) numHaploidSamples;
			haploidFrequencies[1] = 1.0 - haploidFrequencies[0];
		}

		//check if we stop
		maxF = fabs(diploidFrequencies[0] - diploidFrequencies_old[0]);
		double tmp = fabs(diploidFrequencies[1] - diploidFrequencies_old[1]);
		if(tmp > maxF) maxF = tmp;
		tmp = fabs(diploidFrequencies[2] - diploidFrequencies_old[2]);
		if(tmp > maxF) maxF = tmp;
		tmp = fabs(haploidFrequencies[0] - haploidFrequencies_old[0]);
		if(tmp > maxF) maxF = tmp;
		tmp = fabs(haploidFrequencies[1] - haploidFrequencies_old[1]);
		if(tmp > maxF) maxF = tmp;
	}

	alleleFrequency = (numDiploidSamples * (2.0 * diploidFrequencies[2] + diploidFrequencies[1]) + numHaploidSamples * haploidFrequencies[1]) / (2.0 * numDiploidSamples + numHaploidSamples);

	if(alleleFrequency > 0.5) MAF = 1.0 - alleleFrequency;
	else MAF = alleleFrequency;
};

coretools::Log10Probability TGenotypeFrequencies::calculateLog10Likelihood(TPopulationLikehoodLocus & samples){
	return calculateLog10Likelihood(samples.samples(), samples.numSamples());
};

coretools::Log10Probability TGenotypeFrequencies::calculateLog10Likelihood(TSampleLikelihoods* samples, const uint32_t & numSamples){
	double LL = 0.0;
	for(int i = 0; i < numSamples; i++){
		if(!samples[i].isMissing){
			if(samples[i].isHaploid){
				LL += log10((double) (coretools::Probability) samples[i].glfLikelihood_0 * haploidFrequencies[0]
						  + (double) (coretools::Probability) samples[i].glfLikelihood_1 * haploidFrequencies[1]);
			} else {
				LL += log10((double) (coretools::Probability) samples[i].glfLikelihood_0 * diploidFrequencies[0]
						  + (double) (coretools::Probability) samples[i].glfLikelihood_1 * diploidFrequencies[1]
						  + (double) (coretools::Probability) samples[i].glfLikelihood_2 * diploidFrequencies[2]);
			}
		}
	}
	return coretools::Log10Probability(LL);
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
