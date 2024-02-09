/*
 * SFS.cpp
 *
 *  Created on: Apr 10, 2017
 *      Author: wegmannd
 */

#include "SFS.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/algorithms.h"

namespace Simulations {
using coretools::instances::randomGenerator;
using genometools::Base;
//--------------------------------
// Class to store and SFS
//--------------------------------

SFS::SFS(const std::string &filename) {
	coretools::TInputFile in(filename.c_str(), coretools::FileType::NoHeader);

	for (size_t i = 0; i < in.numCols(); ++i) _numChrPerPop.push_back(in.get<size_t>(i));

	//read dimensions
	_numChr = std::accumulate(_numChrPerPop.begin(), _numChrPerPop.end(), 0);

	//add one to each dimension as what is given is #haplotypes = #entries - 1
	_dimensions.resize(_numChrPerPop.size());
	for(size_t i = 0; i < _numChrPerPop.size(); ++i){
		_dimensions[i] = _numChrPerPop[i] + 1;
	}

	//read values
	in.popFront();
	std::vector<double> vec;
	for (size_t i = 0; i < in.numCols(); ++i) vec.push_back(in.get<double>(i));

	// now store as fraction
	coretools::fillFromNormalized(_sfs, vec);
	_sfsPicker.init(_sfs);
}

SFS::SFS(const SFS &other, double MonoFrac) {
	// construct from other SFS, but rescale SFS to have a specific fraction of monomorphic sites
	// now copy and rescale
	_dimensions = other._dimensions;
	_numChrPerPop = other._numChrPerPop;
	_numChr = other._numChr;
	double sum = 0.0;
	_sfs.push_back(MonoFrac);
	for (size_t i = 1; i < other._sfs.size(); ++i) {
		_sfs.push_back(other._sfs[i]);
		sum += other._sfs[i];
	}
	sum = sum / (1.0 - MonoFrac);
	for (auto &s : _sfs) s /= sum;
	_sfsPicker.init(_sfs);
}

SFS::SFS(size_t numChr, double theta) {
	// generate sfs from theta
	_numChrPerPop = { numChr };
	_dimensions = { numChr + 1 };
	_numChr = numChr;
	float sum = 0;
	_sfs.push_back(0);
	for (size_t i = 1; i < _dimensions[0] + 1; ++i) {
		_sfs.push_back(theta / i);
		sum += _sfs.back();
	}
	if (sum > 1.0) UERROR("The choice of theta and sample size results in too many mutations in the SFS!");
	_sfs.front() = 1.0 - sum;
	_sfsPicker.init(_sfs);
}

SFS::SFS(size_t numChr, size_t onlyThisBin) {
	// set all to zero except this one bin
	_numChrPerPop = { numChr };
	_dimensions = { numChr + 1 };
	_numChr = numChr;
	_sfs.resize(numChr + 1, 0.);
	_sfs[onlyThisBin] = 1.0;
	_sfsPicker.init(_sfs);
}

void SFS::writeToFile(const std::string& filename, bool writeLog) const {
	std::ofstream out(filename.c_str());
	if (!out) UERROR("Failed to open file '", filename, "' for writing!");

	//write dimensions
	auto it = _dimensions.begin();
	out << *it;
	for(; it != _dimensions.end(); ++it){
		out << "\t" << *it - 1; //substract one as we write #haplotypes, not #entries
	}
	out << "\n";

	if (writeLog) {
		out << log(_sfs[0]);
		for (size_t i = 1; i < _sfs.size(); ++i) { out << " " << log(_sfs[i]); }
	} else {
		out << _sfs[0];
		for (size_t i = 1; i < _sfs.size(); ++i) { out << " " << _sfs[i]; }
	}
	out << "\n";
	out.close();
}

namespace /* anonymous */ {

constexpr size_t is_odd(size_t x) noexcept { return x & 1; }

}

void SFS::_setDerivedDiploid(size_t l, TSimulatorHaplotypes & haplotypes, size_t N, size_t k, size_t shift, Base derived){
	auto these = _picker.pick(N, k);
	for(auto& i : these){
		size_t index = i + shift;
		size_t ind = index / 2;
		bool hap = is_odd(index);
		haplotypes(ind, hap, l) = derived;
	}
}

void SFS::_setDerivedHaploid(size_t l, TSimulatorHaplotypes & haplotypes, size_t N, size_t k, size_t shift, Base derived){
	//Note: SFS site matches # ind, not # haplotypes. Both haplotypes of an individual are always identical
	auto these = _picker.pick(N, k);
	for(auto& i : these){
		size_t ind = i + shift;
		haplotypes(ind, 0, l) = derived;
		haplotypes(ind, 1, l) = derived;
	}
}

size_t SFS::_simulateSite(size_t l, TSimulatorHaplotypes & haplotypes, Base ancestral, Base derived,
		std::function<void(size_t, TSimulatorHaplotypes &, size_t, size_t, size_t, Base)> func){
		//void (&func)(const size_t, TSimulatorHaplotypes &, const size_t, const size_t, const size_t, const Base)){
	//returns allele count
	//set all as ancestral
	for (size_t i = 0; i < haplotypes.size(); ++i) {
		haplotypes(i, 0, l) = ancestral;
		haplotypes(i, 1, l) = ancestral;
	}

	// pick derived allele frequency
	const auto index = _sfsPicker(randomGenerator().getRand());

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


size_t SFS::simulateSiteHaploid(size_t l, TSimulatorHaplotypes & haplotypes, Base ancestral, Base derived){
	return _simulateSite(l, haplotypes, ancestral, derived,
			[this](const size_t l, TSimulatorHaplotypes & haplotypes, const size_t N, const size_t k, const size_t shift, const Base derived){
				return _setDerivedHaploid(l, haplotypes, N, k, shift, derived);
			});
}

size_t SFS::simulateSiteDiploid(size_t l, TSimulatorHaplotypes & haplotypes, Base ancestral, Base derived){
	return _simulateSite(l, haplotypes, ancestral, derived,
				[this](size_t l, TSimulatorHaplotypes & haplotypes, size_t N, size_t k, size_t shift, Base derived){
					return _setDerivedDiploid(l, haplotypes, N, k, shift, derived);
				});
}

double SFS::calcLLOneSite(const std::vector<double> &gl) {
	double LL = 0.0;
	for (size_t i = 0; i < _sfs.size(); ++i) LL += exp(gl[i]) * _sfs[i];
	return log(LL);
}

} // namespace Simulations
