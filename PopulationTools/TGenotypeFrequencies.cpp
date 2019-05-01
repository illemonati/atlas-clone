/*
 * TGenotypeFreqencies.cpp
 *
 *  Created on: Apr 30, 2019
 *      Author: wegmannd
 */

#include "TGenotypeFrequencies.h"

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

void TGenotypeFrequencies::guess(TPopulationLikehoodStorage & samples){
	guess(samples.samples, samples.numSamples);
};

void TGenotypeFrequencies::guess(TSampleLikelihoods* samples, int numSamples){
	//calculate by using MLE genotype for each individual
	clear();

	for(int i = 0; i < numSamples; i++){
		if(!samples[i].isMissing){
			if(samples[i].isHaploid){
				if(samples[i].phredLikelihood_1 < samples[i].phredLikelihood_0) haploidFrequencies[1] += 1.0;
				else haploidFrequencies[0] += 1.0;
				++numHaploidSamples;
			} else {
				if(samples[i].phredLikelihood_1 < samples[i].phredLikelihood_0){
					if(samples[i].phredLikelihood_2 < samples[i].phredLikelihood_1) diploidFrequencies[2] += 1.0;
					else diploidFrequencies[1] += 1.0;
				} else {
					if(samples[i].phredLikelihood_2 < samples[i].phredLikelihood_0) diploidFrequencies[2] += 1.0;
					else diploidFrequencies[0] += 1.0;
				}
				++numDiploidSamples;
			}
		}
	}

	//normalize
	normalize();
};

void TGenotypeFrequencies::estimate(TPopulationLikehoodStorage & samples, TQualityMap & phredToGTLMap, double epsilonF){
	estimate(samples.samples, samples.numSamples, phredToGTLMap, epsilonF);
};

void TGenotypeFrequencies::estimate(TSampleLikelihoods* samples, int numSamples, TQualityMap & phredToGTLMap, double epsilonF){
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
					weights[0] = phredToGTLMap[samples[i].phredLikelihood_0] * haploidFrequencies_old[0];
					weights[1] = phredToGTLMap[samples[i].phredLikelihood_1] * haploidFrequencies_old[1];

					double sum = weights[0] + weights[1];

					haploidFrequencies[0] += weights[0] / sum;
				} else {
					weights[0] = phredToGTLMap[samples[i].phredLikelihood_0] * diploidFrequencies_old[0];
					weights[1] = phredToGTLMap[samples[i].phredLikelihood_1] * diploidFrequencies_old[1];
					weights[2] = phredToGTLMap[samples[i].phredLikelihood_2] * diploidFrequencies_old[2];

					double sum = weights[0] + weights[1] + weights[2];

					diploidFrequencies[0] += weights[0] / sum;
					diploidFrequencies[2] += weights[2] / sum;
				}
			}
		}

		diploidFrequencies[0] /= (double) numDiploidSamples;
		diploidFrequencies[2] /= (double) numDiploidSamples;
		diploidFrequencies[1] = 1.0 - diploidFrequencies[0] - diploidFrequencies[2];

		haploidFrequencies[0] /= (double) numHaploidSamples;
		haploidFrequencies[1] = 1.0 - haploidFrequencies[0];


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

	//now set allele frequencies
	if(numHaploidSamples > 0)
		alleleFrequency = (numDiploidSamples * (2.0 * diploidFrequencies[2] + diploidFrequencies[1]) + numHaploidSamples * haploidFrequencies[1]) / (2.0 * numDiploidSamples + numHaploidSamples);
	else
		alleleFrequency = (numDiploidSamples * (2.0 * diploidFrequencies[2] + diploidFrequencies[1])) / (2.0 * numDiploidSamples);
	if(alleleFrequency > 0.5) MAF = 1.0 - alleleFrequency;
	else MAF = alleleFrequency;
	std::cout << "MAF " << MAF << std::endl;
};

double TGenotypeFrequencies::calculateLog10Likelihood(TPopulationLikehoodStorage & samples, TQualityMap & phredToGTLMap){
	return calculateLog10Likelihood(samples.samples, samples.numSamples, phredToGTLMap);
};

double TGenotypeFrequencies::calculateLog10Likelihood(TSampleLikelihoods* samples, int numSamples, TQualityMap & phredToGTLMap){
	double LL = 0.0;
	for(int i = 0; i < numSamples; i++){
		if(!samples[i].isMissing){
			if(samples[i].isHaploid){
				LL += log10(phredToGTLMap[samples[i].phredLikelihood_0] * haploidFrequencies[0]
						  + phredToGTLMap[samples[i].phredLikelihood_1] * haploidFrequencies[1]);
			} else {
				LL += log10(phredToGTLMap[samples[i].phredLikelihood_0] * diploidFrequencies[0]
						  + phredToGTLMap[samples[i].phredLikelihood_1] * diploidFrequencies[1]
						  + phredToGTLMap[samples[i].phredLikelihood_2] * diploidFrequencies[2]);
			}
		}
	}
	return LL;
};

