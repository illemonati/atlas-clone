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
#include "debugtools.h"

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
	//void guess(TSampleLikelihoods* samples, int numSamples);
	void estimate(TPopulationLikehoodLocus & samples, const double & epsilonF);
	//void estimate(TSampleLikelihoods* samples, const uint32_t & numSamples, const double & epsilonF);
	coretools::Log10Probability calculateLog10Likelihood(TPopulationLikehoodLocus & samples);
	//coretools::Log10Probability calculateLog10Likelihood(TSampleLikelihoods* samples, const uint32_t & numSamples);
	void writeDiploidFrequencies(coretools::TOutputFile & out);
	void writeHaploidFrequencies(coretools::TOutputFile & out);
	uint32_t numHaploid(){ return numHaploidSamples; };
	uint32_t numDiploid(){ return numDiploidSamples; };

	//template functions to estimate frequencies
	//we use templates so it works with genometools::TSampleLikelihoods (in population tools) and non-standardize likelihoods from GLFs (in major minor)
	template <typename SampleType>
	void guess(SampleType samples, const uint32_t & numSamples){
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

		if(numDiploidSamples > 0){
			diploidFrequencies[0] = (double) diploidCounts[0] / (double) numDiploidSamples;
			diploidFrequencies[1] = (double) diploidCounts[1] / (double) numDiploidSamples;
			diploidFrequencies[2] = (double) diploidCounts[2] / (double) numDiploidSamples;
		} else {
			diploidFrequencies[0] = 0.0; diploidFrequencies[1] = 0.0; diploidFrequencies[2] = 0.0;
		}

		if(numHaploidSamples > 0){
			haploidFrequencies[0] = (double) haploidCounts[0] / (double) numHaploidSamples;
			haploidFrequencies[1] = (double) haploidCounts[1] / (double) numHaploidSamples;
		} else {
			haploidFrequencies[0] = 0.0; haploidFrequencies[1] = 0.0;
		}
	};

	template <typename SampleType>
	void estimate(SampleType samples, const uint32_t & numSamples, const double & epsilonF){
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
				diploidFrequencies_tmp[1] = 1.0 - (diploidFrequencies_tmp[0] + diploidFrequencies_tmp[2]); //1 - sum ensures range despite numeric inaccuracies

				//check if we stop
				maxF = fabs(diploidFrequencies_tmp[0] - diploidFrequencies[0]);
				double tmp = fabs(diploidFrequencies_tmp[1] - diploidFrequencies[1]);
				if(tmp > maxF) maxF = tmp;
				tmp = fabs(diploidFrequencies_tmp[2] - diploidFrequencies[2]);
				if(tmp > maxF) maxF = tmp;

				diploidFrequencies[0] = diploidFrequencies_tmp[0];
				diploidFrequencies[1] = diploidFrequencies_tmp[1];
				diploidFrequencies[2] = diploidFrequencies_tmp[2];
			}

			if(numHaploidSamples > 0){
				haploidFrequencies_tmp[0] /= (double) numHaploidSamples;
				haploidFrequencies_tmp[1] = 1.0 - haploidFrequencies_tmp[0];

				//check if we stop
				double tmp = fabs(haploidFrequencies_tmp[0] - haploidFrequencies[0]);
				if(tmp > maxF) maxF = tmp;
				tmp = fabs(haploidFrequencies_tmp[1] - haploidFrequencies[1]);
				if(tmp > maxF) maxF = tmp;

				haploidFrequencies[0] = haploidFrequencies_tmp[0];
				haploidFrequencies[1] = haploidFrequencies_tmp[1];
			}
		}

		//update parameters
		alleleFrequency = (numDiploidSamples * (2.0 * (double) diploidFrequencies[2] + (double) diploidFrequencies[1]) + numHaploidSamples * (double) haploidFrequencies[1]) / (2.0 * numDiploidSamples + numHaploidSamples);

		if(alleleFrequency > 0.5) MAF = 1.0 - alleleFrequency;
		else MAF = alleleFrequency;
	};

	template <typename SampleType>
	coretools::Log10Probability calculateLog10Likelihood(SampleType samples, const uint32_t & numSamples){
		double LL = 0.0;
		for(int i = 0; i < numSamples; i++){
			if(!samples[i].isMissing()){
				if(samples[i].isHaploid()){
					LL += log10((double) (Probability) samples[i][genometools::haploidFirst] * haploidFrequencies[0]
							  + (double) (Probability) samples[i][genometools::haploidSecond] * haploidFrequencies[1]);
				} else {

					double tmp = (double) (Probability) samples[i][genometools::homoFirst] * diploidFrequencies[0]
							   + (double) (Probability) samples[i][genometools::het] * diploidFrequencies[1]
							   + (double) (Probability) samples[i][genometools::homoSecond] * diploidFrequencies[2];

					LL += log10((double) (Probability) samples[i][genometools::homoFirst] * diploidFrequencies[0]
							  + (double) (Probability) samples[i][genometools::het] * diploidFrequencies[1]
							  + (double) (Probability) samples[i][genometools::homoSecond] * diploidFrequencies[2]);
				}
			}
		}
		return coretools::Log10Probability(LL);
	};

};

}; //end namespace

#endif /* POPULATIONTOOLS_TGENOTYPEFREQUENCIES_H_ */
