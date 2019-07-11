/*
 * TPopulationLikelihoodStorage.h
 *
 *  Created on: Apr 30, 2019
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TPOPULATIONLIKELIHOODLOCUS_H_
#define POPULATIONTOOLS_TPOPULATIONLIKELIHOODLOCUS_H_

#include <stdint.h>
#include <iostream>

//------------------------------------------------
// TSampleLikelihoods
//------------------------------------------------
class TSampleLikelihoods{
public:
	uint16_t glfLikelihood_0;
	uint16_t glfLikelihood_1;
	uint16_t glfLikelihood_2;

	bool isHaploid;
	bool isMissing;

	uint8_t operator[](int genotype){
		if(genotype == 0) return glfLikelihood_0;
		if(genotype == 1) return glfLikelihood_1;
		if(isHaploid) throw "Genotype has to be 0 or 1 for haploid samples!";
		if(genotype == 2) return glfLikelihood_2;
		throw "Genotype has to be 0, 1 or 2 for diploid samples!";
	};
};

//------------------------------------------------
// TPopulationLikehoodStorage
// class used when reading line by line
//------------------------------------------------
class TPopulationLikehoodLocus{
public:
	int numSamples;
	TSampleLikelihoods* samples;

	TPopulationLikehoodLocus(){
		numSamples = 0;
		samples = nullptr;
	};

	TPopulationLikehoodLocus(TPopulationLikehoodLocus & other){
		numSamples = 0;
		resize(other.numSamples);

		//copy data
		for(int i=0; i<numSamples; i++){
			samples[i] = other.samples[i];
		}
	};

	TPopulationLikehoodLocus(int NumSamples){
		numSamples = 0;
		resize(NumSamples);
	};

	~TPopulationLikehoodLocus(){
		clear();
	};

	void clear(){
		if(numSamples > 0){
			delete[] samples;
		}
		numSamples = 0;
	};

	void resize(int NumSamples){
		if(NumSamples != numSamples){
			clear();
			numSamples = NumSamples;
			samples = new TSampleLikelihoods[numSamples];
		}
	};

	TSampleLikelihoods & operator[](int index){
		return samples[index];
	};
};



#endif /* POPULATIONTOOLS_TPOPULATIONLIKELIHOODLOCUS_H_ */
