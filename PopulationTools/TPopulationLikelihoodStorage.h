/*
 * TPopulationLikelihoodStorage.h
 *
 *  Created on: Apr 30, 2019
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TPOPULATIONLIKELIHOODSTORAGE_H_
#define POPULATIONTOOLS_TPOPULATIONLIKELIHOODSTORAGE_H_

#include <stdint.h>

//------------------------------------------------
// TSampleLikelihoods
//------------------------------------------------
class TSampleLikelihoods{
public:
	uint8_t phredLikelihood_0;
	uint8_t phredLikelihood_1;
	uint8_t phredLikelihood_2;

	bool isHaploid;
	bool isMissing;

	uint8_t operator[](int genotype){
		if(genotype == 0) return phredLikelihood_0;
		if(genotype == 1) return phredLikelihood_1;
		if(isHaploid) throw "Genotype has to be 0 or 1 for haploid samples!";
		if(genotype == 2) return phredLikelihood_2;
		throw "Genotype has to be 0, 1 or 2 for diploid samples!";
	};
};

//------------------------------------------------
// TPopulationLikehoodStorage
// class used when reading line by line
//------------------------------------------------
class TPopulationLikehoodStorage{
public:
	int numSamples;
	TSampleLikelihoods* samples;

	TPopulationLikehoodStorage(){
		numSamples = 0;
		samples = nullptr;
	};

	TPopulationLikehoodStorage(int NumSamples){
		numSamples = 0;
		resize(NumSamples);
	};

	~TPopulationLikehoodStorage(){
		clear();
	};

	void clear(){
		if(numSamples > 0)
			delete[] samples;
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



#endif /* POPULATIONTOOLS_TPOPULATIONLIKELIHOODSTORAGE_H_ */
