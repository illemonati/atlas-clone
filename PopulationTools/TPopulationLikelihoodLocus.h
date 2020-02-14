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

	uint16_t operator[](int genotype){
		if(genotype == 0) return glfLikelihood_0;
		if(genotype == 1) return glfLikelihood_1;
		if(isHaploid) throw "Genotype has to be 0 or 1 for haploid samples!";
		if(genotype == 2) return glfLikelihood_2;
		throw "Genotype has to be 0, 1 or 2 for diploid samples!";
	};

	void setAsMissing(){
		glfLikelihood_0 = 0.0;
		glfLikelihood_1 = 0.0;
		glfLikelihood_2 = 0.0;

		isHaploid = false;
		isMissing = true;
	};
};

//------------------------------------------------
// TPopulationLikehoodLocus
// class used when reading line by line
//------------------------------------------------
class TPopulationLikehoodLocus{
private:
	uint32_t _numSamples;
	uint32_t _storageSize;
	TSampleLikelihoods* _samples; //allows simple access to subsets based on populations

	void _copy(const TPopulationLikehoodLocus & other);

public:
	TPopulationLikehoodLocus();
	TPopulationLikehoodLocus(const TPopulationLikehoodLocus & other){
		_copy(other);
	};

	TPopulationLikehoodLocus& operator=(const TPopulationLikehoodLocus & other){
		_copy(other);
		return *this;
	};

	TPopulationLikehoodLocus(int NumSamples);
	~TPopulationLikehoodLocus(){
		clear();
	};

	void clear();
	void resize(uint32_t NumSamples);

	TSampleLikelihoods & operator[](int index){
		return _samples[index];
	};

	TSampleLikelihoods* samples(){
		return _samples;
	};

	uint32_t numSamples(){
		return _numSamples;
	};

	uint32_t numSamplesWithData();
	bool hasData();
	void fillAsMissing();
};

//------------------------------------------------
// TPopulationLikehoodRegion
// class used when reading line by line
//------------------------------------------------
class TPopulationLikehoodWindow{
private:
	uint32_t _numSamples;
	uint32_t _numLoci;
	uint32_t _storageSize;
	TPopulationLikehoodLocus* _loci; //allows simple access to subsets based on populations

	void _init();

public:
	TPopulationLikehoodWindow();
	TPopulationLikehoodWindow(TPopulationLikehoodWindow & other);
	TPopulationLikehoodWindow(uint32_t NumLoci, uint32_t NumSamples);
	TPopulationLikehoodWindow(uint32_t NumLoci);

	~TPopulationLikehoodWindow(){
		clear();
	};

	void clear();
	void resize(uint32_t NumLoci);
	void resize(uint32_t NumLoci, uint32_t NumSamples);

	TPopulationLikehoodLocus & operator[](int locus){
		return _loci[locus];
	};

	uint32_t numSamples(){
		return _numSamples;
	};

	uint32_t numLoci(){
		return _numLoci;
	};

	uint32_t numLociwithData();
	bool hasData();
	void fillAsMissing();
	bool individualHasMissingData(uint32_t individual);
};


#endif /* POPULATIONTOOLS_TPOPULATIONLIKELIHOODLOCUS_H_ */
