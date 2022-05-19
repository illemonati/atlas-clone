/*
 * TAllelicDepthCounts.cpp
 *
 *  Created on: Feb 10, 2020
 *      Author: wegmannd
 */

#include "TAllelicDepthCounts.h"

#include <algorithm>
#include <cstdint>
#include <map>
#include <ostream>
#include <vector>

#include "GenotypeTypes.h"
#include "TFile.h"
#include "TGenotypeData.h"
#include "TParameters.h"
#include "TSite.h"
#include "TWindow.h"
#include "stringFunctions.h"

namespace GenomeTasks{
using coretools::instances::parameters;
using coretools::instances::logfile;

namespace impl {
constexpr size_t index(size_t i1, size_t i2, size_t i3, size_t i4, size_t N) noexcept {
	assert(i1 < N && i2 < N && i3 < N && i4 < N);
	return N * (N * (N * i1 + i2) + i3) + i4;
}
} // namespace impl

TAllelicDepthCounts::TAllelicDepthCounts(){
	_initialized = false;
	_maxAllelicDepth = 0;
	_size = 0;
};

TAllelicDepthCounts::TAllelicDepthCounts(const uint32_t MaxAllelicDepth){
	_initialized = false;
	resize(MaxAllelicDepth);
};

void TAllelicDepthCounts::resize(const uint32_t MaxAllelicDepth){
	if(_maxAllelicDepth != MaxAllelicDepth){
		_freeStorage();
		_maxAllelicDepth = MaxAllelicDepth;
		_size = _maxAllelicDepth + 1;
		_counts.resize(_size*_size*_size*_size);

		//allocate
		/*_counts = new uint32_t***[_size];
		for(uint32_t i=0; i<_size; ++i){
			_counts[i] = new uint32_t**[_size];
			for(uint32_t j=0; j<_size; ++j){
				_counts[i][j] = new uint32_t*[_size];
				for(uint32_t k=0; k<_size; ++k){
					_counts[i][j][k] = new uint32_t[_size];
				}
			}
			}*/
		_initialized = true;
	}

	//set all counts to zero
	clear();
};

void TAllelicDepthCounts::clear(){
	std::fill(_counts.begin(), _counts.end(), 0);
};

void TAllelicDepthCounts::_freeStorage(){
	_counts.clear();
};


void TAllelicDepthCounts::addSite(const GenotypeLikelihoods::TBaseCounts & alleleCounts){
	using genometools::Base;
	const auto aA = alleleCounts[Base::A];
	const auto aC = alleleCounts[Base::C];
	const auto aG = alleleCounts[Base::G];
	const auto aT = alleleCounts[Base::T];
	if (aA < _size && aC < _size && aG < _size && aT < _size) {
		++_counts[impl::index(aA, aC, aG, aT, _size)];
	}
};

void TAllelicDepthCounts::addSiteZeroDepth(){
	++_counts.front();
};

void TAllelicDepthCounts::write(const std::string &filename, bool printEmpty){
	//open file
	coretools::TOutputFile out(filename);

	//write header
	out.writeHeader({"A", "C", "G", "T", "Depth", "majorAllele", "majorDepth","minorAllele", "minorDepth","Counts"});

	//write counts
	//uint32_t max = 0;
	for(uint32_t i=0; i<_size; ++i){
		for(uint32_t j=0; j<_size; ++j){
			for(uint32_t k=0; k<_size; ++k){
				for(uint32_t l=0; l<_size; ++l){
					if(printEmpty || _counts[impl::index(i,j,k,l, _size)] > 0){
						//write numA, C, G and T and depth
						out << i << j << k << l << i+j+k+l;

						//find max
						uint32_t max = i;
						if(j > max) max = j;
						if(k > max) max = k;
						if(l > max) max = l;

						//identify those that are at max
						std::vector<char> tmp;
						if(i == max) tmp.push_back('A');
						if(j == max) tmp.push_back('C');
						if(k == max) tmp.push_back('G');
						if(l == max) tmp.push_back('T');

						//write major
						out << tmp[0] << max;

						//find minor
						if(tmp.size() > 1){
							out << tmp[1] << 0;
						} else {
							//find second
							uint32_t second = 0;
							if(i < max && i > second) second = i;
							if(j < max && j > second) second = j;
							if(k < max && k > second) second = k;
							if(l < max && l > second) second = l;

							//print minor
							if(i == second) out << 'A' << i;
							else if(j == second) out << 'C' << j;
							else if(k == second) out << 'G' << k;
							else out << 'T' << l;
						}

						//write counts
						out << _counts[impl::index(i,j,k,l, _size)] << std::endl;
					}
				}
			}
		}
	}
};

//------------------------------------------
// TAllelicDepth
//------------------------------------------
TAllelicDepth::TAllelicDepth() : TGenome_windows(){
	logfile().list("Will assemble allelic depth up to a max depth of " + coretools::str::toString(_readUpToDepth) + ". (parameter 'maxDepth')");
	if(_readUpToDepth > 100){
		logfile().warning("Allocating count table for a max depth of " + coretools::str::toString(_readUpToDepth) + " uses a lot of memory! Use argument maxDepth to limit.");
	}
	_counts.resize(_readUpToDepth);

	if(parameters().parameterExists("printAll")){
		_writeEmpty = true;
		logfile().list("Will write full table, including cells with zero counts. (parameter 'printAll')");
	} else {
		_writeEmpty = false;
		logfile().list("Will only print cells with non-zero counts. (use 'printAll' to print all cells)");
	}
};

void TAllelicDepth::_handleWindow(){
	logfile().listFlushTime("Adding sites to allelic depth table ...");
	GenotypeLikelihoods::TBaseCounts alleleCounts;
	for(auto& s : _window){
		s.countAlleles(alleleCounts);
		_counts.addSite(alleleCounts);
	}
	logfile().doneTime();
};

void TAllelicDepth::quantifyAlleleicDepth(){
	_traverseBAMWindows();

	//write to file
	std::string outputFileName = _outputName + "_allelicDepth.txt.gz";
	logfile().listFlush("Writing allelic depth table to '" + outputFileName + "' ...");
	_counts.write(outputFileName, _writeEmpty);
	logfile().done();
};

}; // end namespace
