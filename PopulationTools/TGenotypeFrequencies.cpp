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

void TGenotypeFrequencies::guess(TSampleLikelihoods* samples, int numSamples){
	//calculate by using MLE genotype for each individual
	clear();

	std::array<uint32_t, 3> diploidCounts{};
	std::array<uint32_t, 2> haploidCounts{};

	for(int i = 0; i < numSamples; i++){
		if(!samples[i].isMissing()){
			if(samples[i].isHaploid()){
				if(samples[i][genometools::haploidSecond] < samples[i][genometools::haploidFirst]) ++haploidCounts[1];
				else ++haploidCounts[0];
			} else {
				if(samples[i][genometools::het] < samples[i][genometools::homoFirst]){
					if(samples[i][genometools::homoSecond] < samples[i][genometools::het]) ++diploidCounts[2];
					else ++diploidCounts[1];
				} else {
					if(samples[i][genometools::homoSecond] < samples[i][genometools::homoFirst]) ++diploidCounts[2];
					else ++diploidCounts[0];
				}
			}
		}
	}

	numDiploidSamples = diploidCounts[0] +diploidCounts[1] + diploidCounts[2];
	numHaploidSamples = haploidCounts[0] + haploidCounts[1];

	diploidFrequencies[0] = (double) diploidCounts[0] / (double) numDiploidSamples;
	diploidFrequencies[1] = (double) diploidCounts[1] / (double) numDiploidSamples;
	diploidFrequencies[2] = (double) diploidCounts[2] / (double) numDiploidSamples;
	haploidFrequencies[0] = (double) haploidCounts[0] / (double) numHaploidSamples;
	haploidFrequencies[1] = (double) haploidCounts[1] / (double) numHaploidSamples;
};

void TGenotypeFrequencies::estimate(TPopulationLikehoodLocus & samples, const double  &epsilonF){
	estimate(samples.samples(), samples.numSamples(), epsilonF);
};

void TGenotypeFrequencies::estimate(TSampleLikelihoods* samples, const uint32_t & numSamples, const double & epsilonF){
	//prepare variables
	double weights[3];

	//estimate initial frequencies from MLEs
	guess(samples, numSamples);
	ensureAllFrequenciesAreNonZero();

	//run EM for max 1000 steps
	double maxF = 10000.0;
	int s = 0;

	while(maxF > epsilonF && s < 1000){
		++s;

		//set genofreq
		std::array<double, 3> diploidFrequencies_tmp{};
		std::array<double, 2> haploidFrequencies_tmp{};

		//estimate new genotype frequencies
		for(int i = 0; i < numSamples; i++){
			if(!samples[i].isMissing()){
				if(samples[i].isHaploid()){
					weights[0] = (Probability) samples[i][genometools::haploidFirst] * haploidFrequencies[0];
					weights[1] = (Probability) samples[i][genometools::haploidSecond] * haploidFrequencies[1];

					double sum = weights[0] + weights[1];

					haploidFrequencies_tmp[0] += weights[0] / sum;
				} else {
					weights[0] = (Probability) samples[i][genometools::homoFirst] * diploidFrequencies[0];
					weights[1] = (Probability) samples[i][genometools::het] * diploidFrequencies[1];
					weights[2] = (Probability) samples[i][genometools::homoSecond] * diploidFrequencies[2];

					double sum = weights[0] + weights[1] + weights[2];

					diploidFrequencies_tmp[0] += weights[0] / sum;
					diploidFrequencies_tmp[2] += weights[2] / sum;
				}
			}
		}

		if(numDiploidSamples > 0){
			diploidFrequencies_tmp[0] /= (double) numDiploidSamples;
			diploidFrequencies_tmp[2] /= (double) numDiploidSamples;
			diploidFrequencies_tmp[1] = 1.0 - diploidFrequencies_tmp[0] - diploidFrequencies_tmp[2];
		}

		if(numHaploidSamples > 0){
			haploidFrequencies_tmp[0] /= (double) numHaploidSamples;
			haploidFrequencies_tmp[1] = 1.0 - haploidFrequencies_tmp[0];
		}

		//check if we stop
		maxF = fabs(diploidFrequencies_tmp[0] - diploidFrequencies[0]);
		double tmp = fabs(diploidFrequencies_tmp[1] - diploidFrequencies[1]);
		if(tmp > maxF) maxF = tmp;
		tmp = fabs(diploidFrequencies_tmp[2] - diploidFrequencies[2]);
		if(tmp > maxF) maxF = tmp;
		tmp = fabs(haploidFrequencies_tmp[0] - haploidFrequencies[0]);
		if(tmp > maxF) maxF = tmp;
		tmp = fabs(haploidFrequencies_tmp[1] - haploidFrequencies[1]);
		if(tmp > maxF) maxF = tmp;

		//update parameters
		diploidFrequencies[0] = diploidFrequencies_tmp[0];
		diploidFrequencies[1] = diploidFrequencies_tmp[1];
		diploidFrequencies[2] = diploidFrequencies_tmp[2];
		haploidFrequencies[0] = haploidFrequencies_tmp[0];
		haploidFrequencies[1] = haploidFrequencies_tmp[1];
	}

	//update parameters

	alleleFrequency = (numDiploidSamples * (2.0 * (double) diploidFrequencies[2] + (double) diploidFrequencies[1]) + numHaploidSamples * (double) haploidFrequencies[1]) / (2.0 * numDiploidSamples + numHaploidSamples);

	if(alleleFrequency > 0.5) MAF = 1.0 - alleleFrequency;
	else MAF = alleleFrequency;
};

coretools::Log10Probability TGenotypeFrequencies::calculateLog10Likelihood(TPopulationLikehoodLocus & samples){
	return calculateLog10Likelihood(samples.samples(), samples.numSamples());
};

coretools::Log10Probability TGenotypeFrequencies::calculateLog10Likelihood(TSampleLikelihoods* samples, const uint32_t & numSamples){
	double LL = 0.0;
	for(int i = 0; i < numSamples; i++){
		if(!samples[i].isMissing()){
			if(samples[i].isHaploid()){
				LL += log10((double) (Probability) samples[i][genometools::haploidFirst] * haploidFrequencies[0]
						  + (double) (Probability) samples[i][genometools::haploidSecond] * haploidFrequencies[1]);
			} else {
				LL += log10((double) (Probability) samples[i][genometools::homoFirst] * diploidFrequencies[0]
						  + (double) (Probability) samples[i][genometools::het] * diploidFrequencies[1]
						  + (double) (Probability) samples[i][genometools::homoSecond] * diploidFrequencies[2]);
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
