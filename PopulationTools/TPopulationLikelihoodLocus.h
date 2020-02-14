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

	void _copy(const TPopulationLikehoodLocus & other){
		_numSamples = 0;
		_storageSize = 0;
		resize(other._numSamples);

		//copy data
		for(uint32_t i=0; i<_numSamples; i++){
			_samples[i] = other._samples[i];
		}
	};

public:
	TPopulationLikehoodLocus(){
		_numSamples = 0;
		_storageSize = 0;
		_samples = nullptr;
	};

	TPopulationLikehoodLocus(const TPopulationLikehoodLocus & other){
		_copy(other);
	};

	TPopulationLikehoodLocus& operator=(const TPopulationLikehoodLocus & other){
		_copy(other);
		return *this;
	};

	TPopulationLikehoodLocus(int NumSamples){
		_numSamples = 0;
		_storageSize = 0;
		resize(NumSamples);
	};

	~TPopulationLikehoodLocus(){
		clear();
	};

	void clear(){
		if(_storageSize > 0){
			delete[] _samples;
		}
		_numSamples = 0;
		_storageSize = 0;
	};

	void resize(uint32_t NumSamples){
		if(_storageSize < NumSamples){
			clear();
			_numSamples = NumSamples;
			_storageSize = _numSamples;
			_samples = new TSampleLikelihoods[_storageSize];
		}
	};

	TSampleLikelihoods & operator[](int index){
		return _samples[index];
	};

	TSampleLikelihoods* samples(){
		return _samples;
	};

	uint32_t numSamples(){
		return _numSamples;
	};

	uint32_t numSamplesWithData(){
		uint32_t n = 0;
		for(uint32_t i=0; i<_numSamples; i++){
			if(!_samples[i].isMissing){
				++n;
			}
		}
		return n;
	};

	bool hasData(){
		for(uint32_t i=0; i<_numSamples; i++){
			if(!_samples[i].isMissing){
				return true;
			}
		}
		return false;
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

public:
	TPopulationLikehoodWindow(){
		_numSamples = 0;
		_numLoci = 0;
		_storageSize = 0;
		_loci = nullptr;
	};

	TPopulationLikehoodWindow(TPopulationLikehoodWindow & other){
		_numSamples = 0;
		_storageSize = 0;
		_numLoci = 0;
		resize(other._numSamples);

		//_copy data
		for(uint32_t i=0; i<_numLoci; i++){
			_loci[i] = other._loci[i];
		}
	};

	TPopulationLikehoodWindow(int NumSamples){
		_numSamples = 0;
		_storageSize = 0;
		resize(NumSamples);
	};

	~TPopulationLikehoodWindow(){
		clear();
	};

	void clear(){
		if(_storageSize > 0){
			delete[] _loci;
		}
		_numSamples = 0;
		_storageSize = 0;
	};

	void resize(uint32_t NumLoci){
		if(_storageSize < _numLoci){
			clear();
			_numLoci = NumLoci;
			_storageSize = _numLoci;
			_loci = new TPopulationLikehoodLocus[_storageSize];
		}
	};

	void resize(uint32_t NumLoci, uint32_t NumSamples){
		//resize container
		resize(NumLoci);

		//ensure sample size
		if(NumSamples != _numSamples){
			_numSamples = NumSamples;

			for(uint32_t i=0; i<_numLoci; i++){
				_loci[i].resize(_numSamples);
			}
		}
	};


	TPopulationLikehoodLocus & operator[](int index){
		return _loci[index];
	};

	uint32_t numSamples(){
		return _numSamples;
	};

	uint32_t numLoci(){
		return _numLoci;
	};

	uint32_t numLociwithData(){
		uint32_t n = 0;
		for(uint32_t i=0; i<_numLoci; i++){
			if(_loci[i].hasData()){
				++n;
			}
		}
		return n;
	}

	bool hasData(){
		for(uint32_t i=0; i<_numLoci; i++){
			if(_loci[i].hasData()){
				return true;
			}
		}
		return false;
	};

};


#endif /* POPULATIONTOOLS_TPOPULATIONLIKELIHOODLOCUS_H_ */
