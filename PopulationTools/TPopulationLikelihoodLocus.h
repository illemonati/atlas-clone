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
#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TSampleLikelihoods.h"
#include "TTask.h"

namespace PopulationTools{

using TSampleLikelihoods = genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability>;

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

	TSampleLikelihoods & operator[](const uint32_t & index){
		return _samples[index];
	};

	TSampleLikelihoods* samples(){
		return _samples;
	};

	uint32_t numSamples() const{
		return _numSamples;
	};

	uint32_t numSamplesWithData() const;
	bool hasData() const;
	void fillAsMissingHaploid();
	void fillAsMissingDiploid();

	void print(uint32_t index);
	void print();
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
	void fillAsMissingHaploid();
	void fillAsMissingDiploid();

	bool individualHasMissingData(uint32_t individual);

	void print();
	void print(uint32_t locus);
	void printIndividual(uint32_t ind);
};

}; //end namespace

#endif /* POPULATIONTOOLS_TPOPULATIONLIKELIHOODLOCUS_H_ */
