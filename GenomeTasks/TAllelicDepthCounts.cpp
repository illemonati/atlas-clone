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

#include "genometools/GenotypeTypes.h"
#include "coretools/Files/TFile.h"
#include "TGenotypeData.h"
#include "coretools/Main/TParameters.h"
#include "TSite.h"
#include "TWindow.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenomeTasks{
using coretools::instances::parameters;
using coretools::instances::logfile;

namespace impl {
constexpr size_t index(size_t i1, size_t i2, size_t i3, size_t i4, size_t N) noexcept {
	assert(i1 < N && i2 < N && i3 < N && i4 < N);
	return N * (N * (N * i1 + i2) + i3) + i4;
}
} // namespace impl

TAllelicDepthCounts::TAllelicDepthCounts(size_t MaxAllelicDepth){
	resize(MaxAllelicDepth);
};

void TAllelicDepthCounts::resize(size_t MaxAllelicDepth){
	if(_size != MaxAllelicDepth + 1){
		_size = MaxAllelicDepth + 1;
	}

	//set all counts to zero
	clear();
};

void TAllelicDepthCounts::clear(){
	std::fill(_counts.begin(), _counts.end(), 0);
};

void TAllelicDepthCounts::addSite(const GenotypeLikelihoods::TBaseCounts & alleleCounts){
	using genometools::Base;
	const auto aA = alleleCounts[Base::A];
	const auto aC = alleleCounts[Base::C];
	const auto aG = alleleCounts[Base::G];
	const auto aT = alleleCounts[Base::T];
	if (aA < _size && aC < _size && aG < _size && aT < _size) {
		const auto i = impl::index(aA, aC, aG, aT, _size);
		if (_counts.size() <= i) _counts.resize(i + 1);
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
	for (size_t a = 0; a < _size; ++a) {
		for (size_t c = 0; c < _size; ++c) {
			for (size_t g = 0; g < _size; ++g) {
				for (size_t t = 0; t < _size; ++t) {
					const auto i = impl::index(a, c, g, t, _size);
					if((printEmpty || (_counts.size() > i && _counts[i] > 0))){
						//write numA, C, G and T and depth
						out << a << c << g << t << a+c+g+t;
						//find max
						size_t max = a;
						if(c > max) max = c;
						if(g > max) max = g;
						if(t > max) max = t;

						//identify those that are at max
						std::array<char, 4> tmp;
						size_t size = 0;
						if(a == max) tmp[size++] = 'A';
						if(c == max) tmp[size++] = 'C';
						if(g == max) tmp[size++] = 'G';
						if(t == max) tmp[size++] = 'T';

						//write major
						out << tmp.front() << max;

						//find minor
						if(size > 1){
							out << tmp[1] << 0;
						} else {
							//find second
							size_t second = 0;
							if(a < max && a > second)
								second = a;
							if(c < max && c > second)
								second = c;
							if(g < max && g > second)
								second = g;
							if(t < max && t > second)
								second = t;

							//print minor
							if(a == second) out << 'A' << a;
							else if(c == second) out << 'C' << c;
							else if(g == second) out << 'G' << g;
							else out << 'T' << t;
						}

						//write counts
						out.writeln(_counts[impl::index(a,c,g,t, _size)]);
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
	logfile().list("Will assemble allelic depth up to a max depth of " + coretools::str::toString(_readUpToDepth) + ". (parameter 'readUpToDepth')");
	if(_readUpToDepth > 100){
		logfile().warning("Allocating count table for a max depth of " + coretools::str::toString(_readUpToDepth) + " uses a lot of memory! Use argument readUpToDepth to limit.");
	}

	_counts.resize(_readUpToDepth);

	if(parameters().exists("printAll")){
		_writeEmpty = true;
		logfile().list("Will write full table, including cells with zero counts. (parameter 'printAll')");
	} else {
		_writeEmpty = false;
		logfile().list("Will only print cells with non-zero counts. (use 'printAll' to print all cells)");
	}
};

void TAllelicDepth::_handleWindow(){
	logfile().listFlushTime("Adding sites to allelic depth table ...");
	for(auto& s : _window){
		const auto alleleCounts = s.countAlleles();
		_counts.addSite(alleleCounts);
	}
	logfile().doneTime();
};

void TAllelicDepth::run(){
	_traverseBAMWindows();

	//write to file
	std::string outputFileName = _outputName + "_allelicDepth.txt.gz";
	logfile().listFlush("Writing allelic depth table to '" + outputFileName + "' ...");
	_counts.write(outputFileName, _writeEmpty);
	logfile().done();
};

}; // end namespace
