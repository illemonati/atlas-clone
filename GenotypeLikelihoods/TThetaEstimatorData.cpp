/*
 * TThetaEstimatorData.cpp
 *
 *  Created on: Jun 17, 2018
 *      Author: phaentu
 */

#include "TThetaEstimatorData.h"

#include "TBamFilter.h"
#include "coretools/Containers/TView.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "genometools/GenotypeTypes.h"
#include "coretools/Files/TOutputFile.h"
#include "TGenotypeData.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TSite.h"
#include <algorithm>
#include <math.h>
#include <ostream>
#include <stdio.h>
#include "coretools/devtools.h"

namespace GenotypeLikelihoods {

//-------------------------------------------------------
// TThetaEstimatorTemporaryFile
//-------------------------------------------------------
TThetaEstimatorTemporaryFile::TThetaEstimatorTemporaryFile() { init(""); };

TThetaEstimatorTemporaryFile::TThetaEstimatorTemporaryFile(std::string Filename) { init(Filename); };

void TThetaEstimatorTemporaryFile::init(std::string Filename) {
	filename   = Filename;
	fp         = NULL;
	sizeOfData = sizeof(double) * 10;

	isOpen           = false;
	isOpenForReading = false;
	isOpenForWriting = false;
	wasWritten       = false;
};

void TThetaEstimatorTemporaryFile::openForWriting() {
	if (sizeOfData == 0) UERROR("Can not open temporary data file for theta: file was not initialized!");

	// if file was written, remove it
	clean();

	// now open
	fp               = gzopen(filename.c_str(), "wb");
	isOpenForWriting = true;
	isOpenForReading = false;
	wasWritten       = true;

	isOpen = true;
};

void TThetaEstimatorTemporaryFile::openForReading() {
	if (!wasWritten) UERROR("Can not parse temporary file: file was never written!");

	// make sure file is closed
	close();

	// now open
	fp               = gzopen(filename.c_str(), "rb");
	isOpenForWriting = false;
	isOpenForReading = true;

	isOpen = true;
};

void TThetaEstimatorTemporaryFile::close() {
	if (isOpen) {
		gzclose(fp);
		isOpen           = false;
		isOpenForWriting = false;
		isOpenForReading = false;
	}
};

void TThetaEstimatorTemporaryFile::clean() {
	close();
	if (wasWritten) {
		remove(filename.c_str());
		wasWritten = false;
	}
};

bool TThetaEstimatorTemporaryFile::isEOF() {
	if (!isOpenForReading) return true;
	return gzeof(fp);
}

void TThetaEstimatorTemporaryFile::save(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) {
	if (!isOpenForWriting) UERROR("Can not add data to '", filename, "': file is closed!");

	gzwrite(fp, genoLik.data(), sizeOfData);
};

bool TThetaEstimatorTemporaryFile::read(GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) {
	if (!isOpenForReading) UERROR("Can not read data from '", filename, "': file is closed!");
	if (gzread(fp, genoLik.data(), sizeOfData) != sizeOfData) {
		// is end-of-file?
		if (gzeof(fp)) return false;

		// is error
		UERROR("Failed to read data from temporary file!");
	}
	return true;
};

//-------------------------------------------------------
// TThetaEstimatorData
//-------------------------------------------------------
TThetaEstimatorData::TThetaEstimatorData() {
	numSitesCoveredTwiceOrMore = 0;
	totNumSitesAdded           = 0;
	numSitesWithData           = 0;
	cumulativeDepth            = 0.0;
	tmpBaseFreq.fill(0.0);
	readState      = false;
	curSite        = 0;
	curRep         = 0;
};

void TThetaEstimatorData::clear() {
	tmpBaseFreq.fill(0.0);
	numSitesCoveredTwiceOrMore = 0;
	totNumSitesAdded           = 0;
	numSitesWithData           = 0;
	cumulativeDepth            = 0.0;
	emptyStorage();
	clearBootstrap();
};

void TThetaEstimatorData::add(const GenotypeLikelihoods::TSite &site,
							  const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) {
	// assumes that emission probabilities were calculated!!
	++totNumSitesAdded;

	// add if site has data
	if (!site.empty()) {
		++numSitesWithData;
		cumulativeDepth += site.depth();

		saveSite(genoLik);

		// add site to base frequency estimation
		site.addToBaseFrequencies(tmpBaseFreq);

		// count sites covered >=2
		if (site.depth() > 1) ++numSitesCoveredTwiceOrMore;
	}
};

TBaseProbabilities TThetaEstimatorData::baseFrequencies() {
	// estimate base frequencies
	return TBaseProbabilities::normalize(tmpBaseFreq);
};

GenotypeLikelihoods::TGenotypeData TThetaEstimatorData::P_G(const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype) {
	using genometools::Genotype;
	// assumes that pGenotype is set!
	GenotypeLikelihoods::TGenotypeData P_G(0.);

	// calculate P_g for each site
	begin();
	do {
		const auto P_g_oneSite = posterior(curGenotypeLikelihoods(), pGenotype);
		std::transform(P_G.begin(), P_G.end(), P_g_oneSite.begin(), P_G.begin(), std::plus<>());

	} while (next());
	return P_G;
};

double TThetaEstimatorData::calcLogLikelihood(const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype) {
	double LL = 0.0;
	begin();
	do {
		LL += log(weightedSum(curGenotypeLikelihoods(), pGenotype));
	} while (next());

	return LL;
};

void TThetaEstimatorData::addToHeader(std::vector<std::string> &header, const std::string prefix) {
	header.push_back(prefix + "depth");
	header.push_back(prefix + "fracMissing");
	header.push_back(prefix + "fracTwoOrMore");
};

void TThetaEstimatorData::writeSite(coretools::TOutputFile &out) {
	if (_isBootstrapped()) {
		out << "NA";
		out << (double)(totNumSitesAdded - numSitesWithData) / (double)totNumSitesAdded;
		out << "NA";
		// out << "NA"; //TODO: check if this NA is needed.
	} else {
		out << cumulativeDepth / (double)totNumSitesAdded;
		out << (double)(totNumSitesAdded - numSitesWithData) / (double)totNumSitesAdded;
		out << (double)numSitesCoveredTwiceOrMore / (double)totNumSitesAdded;
	}
};

void TThetaEstimatorData::bootstrap() {
	// make sure we start empty
	numBootstrapRepsPerEntry.assign(numSitesWithData, 0);

	for (size_t _ = 0; _ < numSitesWithData; ++_) {
		++numBootstrapRepsPerEntry[coretools::instances::randomGenerator().sample(numSitesWithData)];
	}
};

void TThetaEstimatorData::clearBootstrap() {
	numBootstrapRepsPerEntry.clear();
};

bool TThetaEstimatorData::begin() {
	curSite = 0; // first site is at index zero
	curRep  = 1; // index starts at one

	_begin();

	if (_isBootstrapped()) {
		while (readState && numBootstrapRepsPerEntry[curSite] == 0) { readNext(); }
	}

	return (readState);
};

bool TThetaEstimatorData::next() {
	if (!_isBootstrapped())
		readNext();
	else {
		if (curRep < numBootstrapRepsPerEntry[curSite]) {
			++curRep;
			return true;
		} else {
			curRep = 1;
			readNext();
			while (readState && numBootstrapRepsPerEntry[curSite] == 0) { readNext(); }
		}
	}

	return (readState);
};

//-------------------------------------------------------
// TThetaEstimatorDataVector
//-------------------------------------------------------
TThetaEstimatorDataVector::TThetaEstimatorDataVector() : TThetaEstimatorData(){};

void TThetaEstimatorDataVector::saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) {
	// store emission probabilities
	sites.push_back(genoLik);
};

void TThetaEstimatorDataVector::emptyStorage() { sites.clear(); };

void TThetaEstimatorDataVector::readNext() {
	++siteIt;
	++curSite;

	// readState = (siteIt == sites.end())

	if (siteIt == sites.end()) readState = false;
}

void TThetaEstimatorDataVector::_begin() {
	siteIt    = sites.begin();
	readState = true;
};

GenotypeLikelihoods::TGenotypeLikelihoods &TThetaEstimatorDataVector::curGenotypeLikelihoods() { return *siteIt; }

//-------------------------------------------------------
// TThetaEstimatorDataFile
//-------------------------------------------------------
TThetaEstimatorDataFile::TThetaEstimatorDataFile(std::string TmpFileName) : TThetaEstimatorData() {
	dataFileName = TmpFileName;
	sites.init(dataFileName);
	sites.openForWriting();
};

void TThetaEstimatorDataFile::emptyStorage() { sites.clean(); };

void TThetaEstimatorDataFile::saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) {
	sites.save(genoLik);
};

void TThetaEstimatorDataFile::readNext() {
	++curSite;
	if (curSite >= numSitesWithData)
		readState = false;
	else
		readState = sites.read(genotypeLikelihoods);
}

void TThetaEstimatorDataFile::_begin() {
	sites.openForReading();
	--curSite;
	readNext(); // read first! This is required to match begin() of a vector
}

GenotypeLikelihoods::TGenotypeLikelihoods &TThetaEstimatorDataFile::curGenotypeLikelihoods() {
	return genotypeLikelihoods;
};

}; // namespace GenotypeLikelihoods
