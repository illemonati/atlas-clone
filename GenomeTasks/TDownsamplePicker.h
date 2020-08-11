/*
 * TDownsamplePicker.h
 *
 *  Created on: Aug 10, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TDOWNSAMPLEPICKER_H_
#define GENOMETASKS_TDOWNSAMPLEPICKER_H_

#include "TRandomGenerator.h"
#include <set>


//---------------------------------------
// TDownsampleSolutions
//---------------------------------------
class TDownsampleCase{
private:
	uint8_t _N, _k;

public:
	TDownsampleCase(const uint8_t & N, const uint8_t & k){
		_N = N;
		_k = k;
	};

	uint8_t N() const { return _N; };
	uint8_t k() const { return _k; };

	bool operator<(const TDownsampleCase & other) const{
		if(_N < other._N){
			return true;
		} else if(_N > other._N){
			return false;
		} else if(_k < other._k){
			return true;
		} else {
			return false;
		}
	}
};

//---------------------------------------
// TDownsampleSolutions
//---------------------------------------
class TDownsampleSolutions{
private:
	TDownsampleCase _downsampleCase;
	std::vector< std::vector<uint8_t> > _solutions;

public:
	TDownsampleSolutions(const TDownsampleCase & downsamplingCase);

	bool operator<(const TDownsampleSolutions & other) const;
	const std::vector<uint8_t> pick(TRandomGenerator & RandomGenerator) const;

};

//---------------------------------------
// TBionomCoefStorage
//---------------------------------------
class TBionomCoefStorageOneN{
private:
	std::vector<double> _logBinomCoeff;

public:
	TBionomCoefStorageOneN(const uint32_t & N);
	double operator[](const uint32_t & k) const;
};

class TBionomCoefStorage{
private:
	uint32_t _maxNStored_plusOne;
	std::vector<TBionomCoefStorageOneN> _storage;
public:
	TBionomCoefStorage(const uint32_t & MaxN=500);
	double operator[](const TDownsampleCase & dc) const;
};

//---------------------------------------
// TDownsamplePicker
//---------------------------------------
class TDownsamplePicker{
private:
	TRandomGenerator* _randomGenerator;
	uint8_t _maxNStored;

	//storage for precalculated solutions
	mutable std::set<TDownsampleSolutions, std::less<> > _solutions;

	//variables for random picking
	mutable std::vector<uint8_t> _tmpSolution; //used when N > _maxNStored
	TBionomCoefStorage _binomCoefStorage;



public:
	TDownsamplePicker(TRandomGenerator* RandomGenerator, const uint8_t & MaxStored=20);

	const std::vector<uint8_t> pick(const uint32_t & N, const uint32_t & k) const;
};



#endif /* GENOMETASKS_TDOWNSAMPLEPICKER_H_ */
