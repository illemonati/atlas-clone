/*
 * TDownsamplePicker.cpp
 *
 *  Created on: Aug 10, 2020
 *      Author: phaentu
 */

#include "TDownsamplePicker.h"

//---------------------------------------
// TDownsampleSolutions
//---------------------------------------
TDownsampleSolutions::TDownsampleSolutions(const TDownsampleCase & DownsamplingCase):_downsampleCase(DownsamplingCase){
	//fill solutions
	uint32_t numSolutions = choose(_downsampleCase.N(), _downsampleCase.k());
	_solutions.resize(numSolutions);

	//TODO!!
};

bool TDownsampleSolutions::operator<(const TDownsampleSolutions & other) const{
	return _downsampleCase < other._downsampleCase;
};

const std::vector<uint8_t> TDownsampleSolutions::pick(TRandomGenerator & RandomGenerator) const{
	return _solutions[RandomGenerator.pickOne(_solutions.size())];
};


//---------------------------------------
// TBionomCoefStorage
//---------------------------------------
TBionomCoefStorageOneN::TBionomCoefStorageOneN(const uint32_t & N){
	_logBinomCoeff.resize(N+1);
	for(uint32_t k = 0; k<=N; ++k){
		_logBinomCoeff[k] = chooseLog(N, k);
	}
};

double TBionomCoefStorageOneN::operator[](const uint32_t & k) const{
	return _logBinomCoeff[k];
};

TBionomCoefStorage::TBionomCoefStorage(const uint32_t & MaxN){
	_maxNStored_plusOne = MaxN + 1;

	_storage.reserve(_maxNStored_plusOne);
	for(uint32_t N = 0; N<_maxNStored_plusOne; ++N){
		_storage.emplace_back(N);
	}

};

double TBionomCoefStorage::operator[](const TDownsampleCase & dc) const{
	if(dc.N() < _maxNStored_plusOne){
		return _storage[dc.N()][dc.k()];
	} else {
		return chooseLog(dc.N(), dc.k());
	}
};

//---------------------------------------
// TDownsamplePicker
//---------------------------------------
TDownsamplePicker::TDownsamplePicker(TRandomGenerator* RandomGenerator, const uint8_t & MaxStored){
	_randomGenerator = RandomGenerator;
	_maxNStored = MaxStored;
};

const std::vector<uint8_t> TDownsamplePicker::pick(const uint32_t & N, const uint32_t & k) const{
	if(N < _maxNStored){
		//check if case has been calculated before
		TDownsampleCase dc(N, k);
		const auto c = _solutions.emplace(dc).first;
		return c->pick(*_randomGenerator);
	} else {
		//get random combination

	}
};



