/*
 * TPopulationLikelihoodStorage.cpp
 *
 *  Created on: Apr 30, 2019
 *      Author: wegmannd
 */



/*
 * TPopulationLikelihoodStorage.h
 *
 *  Created on: Apr 30, 2019
 *      Author: wegmannd
 */


#include <stdint.h>
#include <iostream>
#include "TPopulationLikelihoodLocus.h"

namespace PopulationTools{

//------------------------------------------------
// TPopulationLikehoodLocus
// class used when reading line by line
//------------------------------------------------
void TPopulationLikehoodLocus::_copy(const TPopulationLikehoodLocus & other){
	_numSamples = 0;
	_storageSize = 0;
	resize(other._numSamples);

	//copy data
	for(uint32_t i=0; i<_numSamples; i++){
		_samples[i] = other._samples[i];
	}
};

TPopulationLikehoodLocus::TPopulationLikehoodLocus(){
	_numSamples = 0;
	_storageSize = 0;
	_samples = nullptr;
};


TPopulationLikehoodLocus::TPopulationLikehoodLocus(int NumSamples){
	_numSamples = 0;
	_storageSize = 0;
	resize(NumSamples);
};

void TPopulationLikehoodLocus::clear(){
	if(_storageSize > 0){
		delete[] _samples;
	}
	_numSamples = 0;
	_storageSize = 0;
};

void TPopulationLikehoodLocus::resize(uint32_t NumSamples){
	if(_storageSize < NumSamples){
		clear();
		_numSamples = NumSamples;
		_storageSize = _numSamples;
		_samples = new TSampleLikelihoods[_storageSize];
	}
};

uint32_t TPopulationLikehoodLocus::numSamplesWithData() const{
	uint32_t n = 0;
	for(uint32_t i=0; i<_numSamples; i++){
		if(!_samples[i].isMissing){
			++n;
		}
	}
	return n;
};

bool TPopulationLikehoodLocus::hasData() const{
	for(uint32_t i=0; i<_numSamples; i++){
		if(!_samples[i].isMissing){
			return true;
		}
	}
	return false;
};

void TPopulationLikehoodLocus::fillAsMissing(){
	for(uint32_t i=0; i<_numSamples; i++){
		_samples[i].setAsMissing();
	}
};

//------------------------------------------------
// TPopulationLikehoodWindow
// class used when reading line by line
//------------------------------------------------
TPopulationLikehoodWindow::TPopulationLikehoodWindow(){
	_init();
};

TPopulationLikehoodWindow::TPopulationLikehoodWindow(TPopulationLikehoodWindow & other){
	_init();
	resize(other._numSamples);

	//_copy data
	for(uint32_t i=0; i<_numLoci; i++){
		_loci[i] = other._loci[i];
	}
};

TPopulationLikehoodWindow::TPopulationLikehoodWindow(uint32_t NumLoci, uint32_t NumSamples){
	_init();
	resize(NumLoci, NumSamples);
};

TPopulationLikehoodWindow::TPopulationLikehoodWindow(uint32_t NumLoci){
	_init();
	resize(NumLoci);
};

void TPopulationLikehoodWindow::_init(){
	_numSamples = 0;
	_numLoci = 0;
	_storageSize = 0;
	_loci = nullptr;
};

void TPopulationLikehoodWindow::clear(){
	if(_storageSize > 0){
		delete[] _loci;
	}
	_numSamples = 0;
	_storageSize = 0;
};

void TPopulationLikehoodWindow::resize(uint32_t NumLoci){
	resize(NumLoci, _numSamples);
};

void TPopulationLikehoodWindow::resize(uint32_t NumLoci, uint32_t NumSamples){
	//resize container
	if(_storageSize < NumLoci){
		clear();
		_numLoci = NumLoci;
		_storageSize = _numLoci;
		_loci = new TPopulationLikehoodLocus[_storageSize];
	}

	//ensure sample size
	_numSamples = NumSamples;
	for(uint32_t i=0; i<_numLoci; i++){
		_loci[i].resize(_numSamples);
	}
};

uint32_t TPopulationLikehoodWindow::numLociwithData(){
	uint32_t n = 0;
	for(uint32_t i=0; i<_numLoci; i++){
		if(_loci[i].hasData()){
			++n;
		}
	}
	return n;
}

bool TPopulationLikehoodWindow::hasData(){
	for(uint32_t i=0; i<_numLoci; i++){
		if(_loci[i].hasData()){
			return true;
		}
	}
	return false;
};

void TPopulationLikehoodWindow::fillAsMissing(){
	for(uint32_t i=0; i<_numLoci; i++){
		_loci[i].fillAsMissing();
	}
};

bool TPopulationLikehoodWindow::individualHasMissingData(uint32_t individual){
	for(uint32_t i=0; i<_numLoci; i++){
		if(_loci[i][individual].isMissing)
			return true;
	}
	return false;
};

}; //end namespace
