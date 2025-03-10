/*
 * SFS.cpp
 *
 *  Created on: Apr 10, 2017
 *      Author: wegmannd
 */

#include "SFS.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/algorithms.h"
#include "coretools/Strings/stringManipulations.h"
#include "coretools/Strings/fillContainer.h"
#include "genometools/Genotypes/TwoBases.h"

namespace Simulations {
using coretools::instances::randomGenerator;
using genometools::Base;
//--------------------------------
// Class to store and SFS
//--------------------------------

SFS::SFS(const std::string &filename) {
	coretools::TInputFile in(filename.c_str(), coretools::FileType::NoHeader, coretools::str::whitespaces, "//");

	// first line reads "SHAPE=<x,y,z>", where x, y, z are dimensions per pop	
	auto tmp = coretools::str::readBefore(coretools::str::readAfter(in.get(0), "<"), ">");
	coretools::str::fillContainerFromString(tmp, _dimensions, "/");

	//numChr = dimensions - 1: what is given is #entries = #haplotypes + 1 
	_numChrPerPop.resize(_dimensions.size());
	_numChr = 0;
	size_t numCells = 1;
	for(size_t i = 0; i < _dimensions.size(); ++i){
		_numChrPerPop[i] = _dimensions[i] - 1;
		_numChr += _numChrPerPop[i];
		numCells *= _dimensions[i];
	}

	// read linbe with values and check dimensionality
	in.popFront();
	if(in.numCols() != numCells){
		UERROR("Error reading SFS from file '", filename, "': number of entries (", in.numCols(), ") does not match number of entries expected from dimensions (", numCells, ")!");
	}

	//read values
	std::vector<double> vec;
	for (size_t i = 0; i < in.numCols(); ++i){
		vec.push_back(in.get<double>(i));
	} 

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
	for (size_t i = 1; i < _dimensions[0]; ++i) {
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

void SFS::writeToFile(std::string_view Filename, bool WriteLog) const {
	coretools::TOutputFile out(Filename);

	out.writeln("#SHAPE=<" + coretools::str::concatenateString(_dimensions, "/") + ">");

	if (WriteLog) {		
		for (size_t i = 0; i < _sfs.size(); ++i){
			 out.write(log(_sfs[i]));
		}
	} else {
		for (size_t i = 0; i < _sfs.size(); ++i){
			 out.write(_sfs[i]);
		}
	}
	out.endln();
}

namespace impl {

constexpr size_t is_odd(size_t N) noexcept { return N & 1; }

}

void SFS::_setDerivedDiploid(size_t l, TSimulatorHaplotypes & haplotypes, size_t N, size_t k, size_t shift, Base derived){
	auto these = _picker.pick(N, k);
	for(auto& i : these){
		const auto index = i + shift;
		const auto ind = index / 2;
		const auto b1 = first(haplotypes[ind][l]);
		const auto b2 = second(haplotypes[ind][l]);
		if (impl::is_odd(index)) {
			haplotypes[ind][l] = genometools::twoBase(b1, derived);
		} else {
			haplotypes[ind][l] = genometools::twoBase(derived, b2);
		}
	}
}

void SFS::_setDerivedHaploid(size_t l, TSimulatorHaplotypes & haplotypes, size_t N, size_t k, size_t shift, Base derived){
	//Note: SFS site matches # ind, not # haplotypes. Both haplotypes of an individual are always identical
	auto these = _picker.pick(N, k);
	for(auto& i : these){
		size_t ind = i + shift;
		haplotypes[ind][l] = genometools::twoBase(derived, derived);
	}
}

size_t SFS::_simulateSite(size_t l, TSimulatorHaplotypes &haplotypes, Base ancestral, Base derived, bool haplo) {
	// void (&func)(const size_t, TSimulatorHaplotypes &, const size_t, const size_t, const size_t, const Base)){
	//returns allele count
	//set all as ancestral
	for (size_t i = 0; i < haplotypes.size(); ++i) {
		haplotypes[i][l] = genometools::twoBase(ancestral, ancestral);
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
			haplo ? _setDerivedHaploid(l, haplotypes, _numChrPerPop[pop], subscripts[pop], shift, derived)
				  : _setDerivedDiploid(l, haplotypes, _numChrPerPop[pop], subscripts[pop], shift, derived);
		}
		shift += _numChrPerPop[pop];
	}
	return std::accumulate(subscripts.begin(), subscripts.end(), 0);
}

size_t SFS::simulateSiteHaploid(size_t l, TSimulatorHaplotypes & haplotypes, Base ancestral, Base derived){
	return _simulateSite(l, haplotypes, ancestral, derived, true);
}

size_t SFS::simulateSiteDiploid(size_t l, TSimulatorHaplotypes & haplotypes, Base ancestral, Base derived){
	return _simulateSite(l, haplotypes, ancestral, derived, false);
}

double SFS::calcLLOneSite(const std::vector<double> &gl) {
	double LL = 0.0;
	for (size_t i = 0; i < _sfs.size(); ++i) LL += exp(gl[i]) * _sfs[i];
	return log(LL);
}

} // namespace Simulations
