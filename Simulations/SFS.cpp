/*
 * SFS.cpp
 *
 *  Created on: Apr 10, 2017
 *      Author: wegmannd
 */

#include "SFS.h"
#include "TRandomGenerator.h"

namespace Simulations {
using coretools::instances::randomGenerator;
//--------------------------------
// Class to store and SFS
//--------------------------------

SFS::SFS(const std::string &filename) {
	std::ifstream file(filename.c_str());
	if (!file) throw "Failed to open SFS file '" + filename + "'!";
	std::vector<double> vec;
	coretools::str::fillContainerFromLineWhiteSpace(file, vec, true, true);

	// now store as fraction
	float sum = 0;
	for (size_t i = 0; i < vec.size(); ++i) {
		sfs.push_back(vec[i]);
		sum += vec[i];
	}
	for (auto &s : sfs) s /= sum;
	_fillFrequencies();
	_fillCumulative();
}

SFS::SFS(const SFS &other, float MonoFrac) {
	// construct from other SFS, but rescale SFS to have a specific fraction of monomorphic sites
	// now copy and rescale
	float sum = 0.0;
	sfs.push_back(MonoFrac);
	for (uint32_t i = 1; i < other.sfs.size(); ++i) {
		sfs.push_back(other.sfs[i]);
		sum += other.sfs[i];
	}
	sum = sum / (1.0 - MonoFrac);
	for (auto &s : sfs) s /= sum;
	_fillFrequencies();
	_fillCumulative();
}

SFS::SFS(uint32_t numChr, float theta) {
	// generate sfs from theta
	float sum = 0;
	for (uint32_t i = 1; i < numChr + 1; ++i) {
		sfs.push_back(theta / i);
		sum += sfs.back();
	}
	if (sum > 1.0) throw "The choice of theta and sample size results in too many mutations in the SFS!";

	sfs.front() = 1.0 - sum;
	_fillFrequencies();
	_fillCumulative();
}

SFS::SFS(uint32_t numChr, uint32_t onlyThisBin) {
	// set all to zero except this one bin
	sfs.resize(numChr + 1, 0.);
	sfs[onlyThisBin] = 1.0;
	_fillFrequencies();
	_fillCumulative();
}

void SFS::_fillFrequencies() {
	sfsFrequencies[0] = 0.0;
	for (uint32_t i = 1; i < sfs.size(); ++i) sfsFrequencies.push_back((float)i / numChromosomes());
}

void SFS::_fillCumulative() {
	sfsCumulative[0] = sfs[0];
	for (uint32_t i = 1; i < sfs.size(); ++i) sfsCumulative.push_back(sfsCumulative[i - 1] + sfs[i]);
	sfsCumulative.back() = 1.0;
}

void SFS::writeToFile(const std::string &filename, const bool &writeLog) {
	std::ofstream out(filename.c_str());
	if (!out) throw "Failed to open file '" + filename + "' for writing!";

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

double SFS::calcLLOneSite(const std::vector<float> &gl) {
	double LL = 0.0;
	for (uint32_t i = 0; i < sfs.size(); ++i) LL += exp(gl[i]) * sfs[i];
	return log(LL);
}

double SFS::getRandomFrequency() {
	return sfsFrequencies[randomGenerator().pickOne(sfsCumulative)];
}

uint32_t SFS::getRandomAlleleCount() {
	return randomGenerator().pickOne(sfsCumulative);
}

//--------------------------------------
// SFSfolded
//--------------------------------------

double SFSfolded::calcLLOneSite(const std::vector<float> &gl) {
	const auto dimension = sfs.size();
	double LL            = 0.0;
	for (uint32_t i = 0; i < dimension; ++i) LL += exp(gl[i]) * sfs[i];

	int j = dimension - 1;
	for (uint32_t i = dimension; i < 2 * dimension - 1; ++i, --j) LL += exp(gl[i]) * sfs[j];

	return log(LL);
}

} // namespace Simulations
