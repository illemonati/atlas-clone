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
#include "Types.h"

namespace PopulationTools{

//------------------------------------------------
// TSampleLikelihoods
//------------------------------------------------
class TSampleLikelihoods{
public:
	BAM::HighPrecisionPhredIntErrorRate glfLikelihood_0;
	BAM::HighPrecisionPhredIntErrorRate glfLikelihood_1;
	BAM::HighPrecisionPhredIntErrorRate glfLikelihood_2;

	bool isHaploid;
	bool isMissing;

	constexpr BAM::HighPrecisionPhredIntErrorRate& operator[](const uint8_t & genotype){
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

	void print(){
		if(isMissing)
			std::cout << "-";
		else {
			if(isHaploid)
				std::cout << glfLikelihood_0 << "," << glfLikelihood_2;
			else
				std::cout << glfLikelihood_0 << "," << glfLikelihood_1 << "," << glfLikelihood_2;
		}
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
	void fillAsMissing();

	void print(uint32_t index){
		_samples[index].print();
	};

	void print(){
		_samples[0].print();
		for(uint32_t i=1; i<_numSamples; ++i){
			std::cout << '\t';
			_samples[i].print();
		}
	};
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

	void print(){
		for(uint32_t l=0; l<_numLoci; ++l){
			std::cout << "Pos" << l << "\t";
			_loci[l].print();
			std::cout << std::endl;
		}
	};

	void print(uint32_t locus){
		_loci[locus].print();
		std::cout << std::endl;
	};

	void printIndividual(uint32_t ind){
		_loci[0].print(ind);
		for(uint32_t l=1; l<_numLoci; ++l){
			std::cout << '\t';
			_loci[l].print(ind);
		}
		std::cout << std::endl;
	};

};

}; //end namespace

#endif /* POPULATIONTOOLS_TPOPULATIONLIKELIHOODLOCUS_H_ */
