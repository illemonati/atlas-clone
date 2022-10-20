/*
 * SFS.cpp
 *
 *  Created on: Apr 10, 2017
 *      Author: wegmannd
 */

#include "SFS.h"

#include <math.h>
#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/algorithms.h"


namespace Simulations {
using coretools::instances::randomGenerator;
//--------------------------------
// Class to store and SFS
//--------------------------------

SFS::SFS(const std::string &filename) {
	std::ifstream file(filename.c_str());
	if (!file) throw "Failed to open SFS file '" + filename + "'!";

	//read dimensions
	coretools::str::fillContainerFromLineWhiteSpace(file, _numChrPerPop, true, true);
	_numChr = std::accumulate(_numChrPerPop.begin(), _numChrPerPop.end(), 0);

	//add one to each dimension as what is given is #haplotypes = #entries - 1
	_dimensions.resize(_numChrPerPop.size());
	for(size_t i = 0; i < _numChrPerPop.size(); ++i){
		_dimensions[i] = _numChrPerPop[i] + 1;
	}

	//read values
	std::vector<double> vec;
	coretools::str::fillContainerFromLineWhiteSpace(file, vec, true, true);

	// now store as fraction
	coretools::fillFromNormalized(sfs, vec);
	coretools::fillCumulative(sfs, sfsCumulative);
}

SFS::SFS(const SFS &other, float MonoFrac) {
	// construct from other SFS, but rescale SFS to have a specific fraction of monomorphic sites
	// now copy and rescale
	_dimensions = other._dimensions;
	_numChrPerPop = other._numChrPerPop;
	_numChr = other._numChr;
	double sum = 0.0;
	sfs.push_back(MonoFrac);
	for (uint32_t i = 1; i < other.sfs.size(); ++i) {
		sfs.push_back(other.sfs[i]);
		sum += other.sfs[i];
	}
	sum = sum / (1.0 - MonoFrac);
	for (auto &s : sfs) s /= sum;
	coretools::fillCumulative(sfs, sfsCumulative);
}

SFS::SFS(uint32_t numChr, float theta) {
	// generate sfs from theta
	_numChrPerPop = { numChr };
	_dimensions = { numChr + 1 };
	_numChr = numChr;
	float sum = 0;
	sfs.push_back(0);
	for (uint32_t i = 1; i < _dimensions[0] + 1; ++i) {
		sfs.push_back(theta / i);
		sum += sfs.back();
	}
	if (sum > 1.0) throw "The choice of theta and sample size results in too many mutations in the SFS!";
	sfs.front() = 1.0 - sum;
	coretools::fillCumulative(sfs, sfsCumulative);
}

SFS::SFS(uint32_t numChr, uint32_t onlyThisBin) {
	// set all to zero except this one bin
	_numChrPerPop = { numChr };
	_dimensions = { numChr + 1 };
	_numChr = numChr;
	sfs.resize(numChr + 1, 0.);
	sfs[onlyThisBin] = 1.0;
	coretools::fillCumulative(sfs, sfsCumulative);
}

void SFS::writeToFile(const std::string &filename, const bool &writeLog) {
	std::ofstream out(filename.c_str());
	if (!out) throw "Failed to open file '" + filename + "' for writing!";

	//write dimensions
	auto it = _dimensions.begin();
	out << *it;
	for(; it != _dimensions.end(); ++it){
		out << "\t" << *it - 1; //substract one as we write #haplotypes, not #entries
	}
	out << "\n";

	if (writeLog) {
		out << log(sfs[0]);
		for (uint32_t i = 1; i < sfs.size(); ++i) { out << " " << log(sfs[i]); }
	} else {
		out << sfs[0];
		for (uint32_t i = 1; i < sfs.size(); ++i) { out << " " << sfs[i]; }
	}
	out << "\n";
	out.close();
}

namespace /* anonymous */ {

constexpr size_t is_odd(size_t x) noexcept { return x & 1; }

}

void SFS::_setDerivedDiploid(const uint32_t l, TSimulatorHaplotypes & haplotypes, const size_t N, const size_t k, const size_t shift, const Base derived){
	auto these = _picker.pick(N, k);
	for(auto& i : these){
		size_t index = i + shift;
		size_t ind = index / 2;
		bool hap = is_odd(index);
		haplotypes(ind, hap, l) = derived;
	}
}

void SFS::_setDerivedHaploid(const uint32_t l, TSimulatorHaplotypes & haplotypes, const size_t N, const size_t k, const size_t shift, const Base derived){
	//Note: SFS site matches # ind, not # haplotypes. Both haplotypes of an individual are always identical
	auto these = _picker.pick(N, k);
	for(auto& i : these){
		size_t ind = i + shift;
		haplotypes(ind, 0, l) = derived;
		haplotypes(ind, 1, l) = derived;
	}
}

size_t SFS::_simulateSite(const uint32_t l, TSimulatorHaplotypes & haplotypes, const Base ancestral, const Base derived,
		std::function<void(const uint32_t, TSimulatorHaplotypes &, const size_t, const size_t, const size_t, const Base)> func){
		//void (&func)(const uint32_t, TSimulatorHaplotypes &, const size_t, const size_t, const size_t, const Base)){
	//returns allele count
	//set all as ancestral
	for (size_t i = 0; i < haplotypes.size(); ++i) {
		haplotypes(i, 0, l) = ancestral;
		haplotypes(i, 1, l) = ancestral;
	}

	// pick derived allele frequency
	size_t index = randomGenerator().pickOne(sfsCumulative);

	//is it polymorphic?
	if(index == 0){
		//monomorphic
		return 0;
	}

	// is polymorphic
	auto subscripts = coretools::getSubscripts(index, _dimensions);
	size_t shift = 0;
	for(size_t pop = 0; pop < _numChrPerPop.size(); ++pop){
		if(subscripts[pop] > 0){
			func(l, haplotypes, _numChrPerPop[pop], subscripts[pop], shift, derived);
		}
		shift += _numChrPerPop[pop];
	}
	return std::accumulate(subscripts.begin(), subscripts.end(), 0);
}


size_t SFS::simulateSiteHaploid(const uint32_t l, TSimulatorHaplotypes & haplotypes, const Base ancestral, const Base derived){
	return _simulateSite(l, haplotypes, ancestral, derived,
			[this](const uint32_t l, TSimulatorHaplotypes & haplotypes, const size_t N, const size_t k, const size_t shift, const Base derived){
				return _setDerivedHaploid(l, haplotypes, N, k, shift, derived);
			});

	/*
	//returns allele count
	//set all as ancestral
	for (size_t i = 0; i < haplotypes.size(); ++i) {
		haplotypes(i, 0, l) = ancestral;
		haplotypes(i, 1, l) = ancestral;
	}

	// pick derived allele frequency
	size_t index = randomGenerator().pickOne(sfsCumulative);

	//is it polymorphic?
	if(index == 0){
		//monomorphic
		return 0;
	}

	// is polymorphic
	auto subscripts = coretools::getSubscripts(index, _dimensions);
	size_t shift = 0;
	for(size_t pop = 0; pop < _dimensions.size(); ++pop){
		_setDerivedDiploid(l, haplotypes, _dimensions[pop], subscripts[pop], shift, derived);
		shift += _dimensions[pop];
	}
	return std::accumulate(subscripts.begin(), subscripts.end(), 0);
	*/
}

size_t SFS::simulateSiteDiploid(const uint32_t l, TSimulatorHaplotypes & haplotypes, const Base ancestral, const Base derived){
	return _simulateSite(l, haplotypes, ancestral, derived,
				[this](const uint32_t l, TSimulatorHaplotypes & haplotypes, const size_t N, const size_t k, const size_t shift, const Base derived){
					return _setDerivedDiploid(l, haplotypes, N, k, shift, derived);
				});
	/*
	//Note: SFS site matches # ind, not # haplotypes. Both haplotypes of an individual are always identical
	//returns allele count
	//set all as ancestral
	for (size_t i = 0; i < haplotypes.size(); ++i) {
		haplotypes(i, 0, l) = ancestral;
		haplotypes(i, 1, l) = ancestral;
	}

	// pick derived allele frequency
	size_t index = randomGenerator().pickOne(sfsCumulative);

	//is it polymorphic?
	if(index == 0){
		//monomorphic
		return 0;
	}

	// is polymorphic
	auto subscripts = coretools::getSubscripts(index, _dimensions);
	size_t shift = 0;
	for(size_t pop = 0; pop < _dimensions.size(); ++pop){
		_setDerivedDiploid(l, haplotypes, _dimensions[pop], subscripts[pop], shift, derived);
		shift += _dimensions[pop];
	}
	return std::accumulate(subscripts.begin(), subscripts.end(), 0);
	*/
}

double SFS::calcLLOneSite(const std::vector<float> &gl) {
	double LL = 0.0;
	for (uint32_t i = 0; i < sfs.size(); ++i) LL += exp(gl[i]) * sfs[i];
	return log(LL);
}

/*
//--------------------------------------
// SFSfolded
//--------------------------------------
SFSfolded::SFSfolded(uint32_t sampleSize, float theta) {
	// generate sfs from theta
	float sum = 0;
	sfs.push_back(0);
	for (uint32_t i = 1; i < sampleSize + 1; ++i) {
		sfs.push_back((theta / i)+(theta/(2*sampleSize-i)));
		sum += sfs.back();
	}
	if (sum > 1.0) throw "The choice of theta and sample size results in too many mutations in the SFS!";
	sfs.front() = 1.0 - sum;
	_fillFrequencies();
	_fillCumulative();
}

double SFSfolded::calcLLOneSite(const std::vector<float> &gl) {
	const auto dimension = sfs.size();
	double LL            = 0.0;
	for (uint32_t i = 0; i < dimension; ++i) LL += exp(gl[i]) * sfs[i];

	int j = dimension - 1;
	for (uint32_t i = dimension; i < 2 * dimension - 1; ++i, --j) LL += exp(gl[i]) * sfs[j];

	return log(LL);
}
*/
} // namespace Simulations
