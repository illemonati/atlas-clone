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

namespace GenotypeLikelihoods {

//-------------------------------------------------------
// TThetaEstimatorTemporaryFile
//-------------------------------------------------------
void TThetaEstimatorTemporaryFile::openForWriting() {
	if (_sizeOfData == 0) UERROR("Can not open temporary data file for theta: file was not initialized!");

	// if file was written, remove it
	clean();

	// now open
	_fp               = gzopen(_filename.c_str(), "wb");
	_isOpenForWriting = true;
	_isOpenForReading = false;
	_wasWritten       = true;
};

void TThetaEstimatorTemporaryFile::openForReading() {
	if (!_wasWritten) UERROR("Can not parse temporary file: file was never written!");

	// make sure file is closed
	close();

	// now open
	_fp               = gzopen(_filename.c_str(), "rb");
	_isOpenForWriting = false;
	_isOpenForReading = true;
};

void TThetaEstimatorTemporaryFile::close() {
	if (_fp) {
		gzclose(_fp);
		_isOpenForWriting = false;
		_isOpenForReading = false;
	}
};

void TThetaEstimatorTemporaryFile::clean() {
	close();
	if (_wasWritten) {
		remove(_filename.c_str());
		_wasWritten = false;
	}
};

bool TThetaEstimatorTemporaryFile::isEOF() {
	if (!_isOpenForReading) return true;
	return gzeof(_fp);
}

void TThetaEstimatorTemporaryFile::save(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) {
	if (!_isOpenForWriting) UERROR("Can not add data to '", _filename, "': file is closed!");

	gzwrite(_fp, genoLik.data(), _sizeOfData);
};

bool TThetaEstimatorTemporaryFile::read(GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) {
	if (!_isOpenForReading) UERROR("Can not read data from '", _filename, "': file is closed!");
	if (gzread(_fp, genoLik.data(), _sizeOfData) != _sizeOfData) {
		// is end-of-file?
		if (gzeof(_fp)) return false;

		// is error
		UERROR("Failed to read data from temporary file!");
	}
	return true;
};

//-------------------------------------------------------
// TThetaEstimatorData
//-------------------------------------------------------
TThetaEstimatorData::TThetaEstimatorData() {
	_numSites2x = 0;
	_numSites           = 0;
	_numSitesData           = 0;
	_cumulativeDepth            = 0.0;
	_readState      = false;
	_curSite        = 0;
	_curRep         = 0;
};

void TThetaEstimatorData::clear() {
	_numSites2x = 0;
	_numSites           = 0;
	_numSitesData           = 0;
	_cumulativeDepth            = 0.0;
	_emptyStorage();
	_numBootstrapRepsPerEntry.clear();
	_baseFreqs.clear();
};

void TThetaEstimatorData::add(const GenotypeLikelihoods::TSite &site,
							  const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) {
	// assumes that emission probabilities were calculated!!
	++_numSites;

	// add if site has data
	if (!site.empty()) {
		++_numSitesData;
		_cumulativeDepth += site.depth();

		_saveSite(genoLik);

		// add site to base frequency estimation
		_baseFreqs.push_back(site.baseFrequencies());

		// count sites covered >=2
		if (site.depth() > 1) ++_numSites2x;
	}
};

TBaseProbabilities TThetaEstimatorData::baseFrequencies() {
	// estimate base frequencies
	TBaseData bd{};
	_begin();
	do {
		bd += _baseFreqs[_curSite];
	} while (_next());

	return TBaseProbabilities::normalize(bd);
};

GenotypeLikelihoods::TGenotypeData TThetaEstimatorData::P_G(const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype) {
	using genometools::Genotype;
	// assumes that pGenotype is set!
	GenotypeLikelihoods::TGenotypeData P_G(0.);

	// calculate P_g for each site
	_begin();
	do {
		const auto P_g_oneSite = posterior(_GL(), pGenotype);
		std::transform(P_G.begin(), P_G.end(), P_g_oneSite.begin(), P_G.begin(), std::plus<>());
	} while (_next());
	return P_G;
};

double TThetaEstimatorData::calcLogLikelihood(const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype) {
	double LL = 0.0;
	_begin();
	do {
		LL += log(weightedSum(_GL(), pGenotype));
	} while (_next());

	return LL;
};

void TThetaEstimatorData::addToHeader(std::vector<std::string> &header, const std::string &prefix) {
	header.push_back(prefix + "depth");
	header.push_back(prefix + "numSites");
	header.push_back(prefix + "numSitesData");
	header.push_back(prefix + "numSites2x+");
	header.push_back(prefix + "fracMissing");
	header.push_back(prefix + "frac2x+");
};

void TThetaEstimatorData::writeSite(coretools::TOutputFile &out) {
	if (_isBootstrapped()) {
		out.write("NA", _numSites, _numSitesData, "NA", double(_numSites - _numSitesData)/_numSites, "NA");
		// out << "NA"; //TODO: check if this NA is needed.
	} else {
		out.write(_cumulativeDepth/_numSites, _numSites, _numSitesData, _numSites2x,
				  double(_numSites - _numSitesData)/_numSites, double(_numSites2x)/_numSites);
	}
};

void TThetaEstimatorData::bootstrap() {
	// make sure we start empty
	_numBootstrapRepsPerEntry.assign(_numSitesData, 0);

	for (size_t _ = 0; _ < _numSitesData; ++_) {
		++_numBootstrapRepsPerEntry[coretools::instances::randomGenerator().sample(_numSitesData)];
	}
};

bool TThetaEstimatorData::_begin() {
	_curSite = 0; // first site is at index zero
	_curRep  = 1; // index starts at one

	_start();

	if (_isBootstrapped()) {
		while (_readState && _numBootstrapRepsPerEntry[_curSite] == 0) { _readNext(); }
	}

	return _readState;
};

bool TThetaEstimatorData::_next() {
	if (!_isBootstrapped())
		_readNext();
	else {
		if (_curRep < _numBootstrapRepsPerEntry[_curSite]) {
			++_curRep;
			return true;
		} else {
			_curRep = 1;
			_readNext();
			while (_readState && _numBootstrapRepsPerEntry[_curSite] == 0) { _readNext(); }
		}
	}

	return _readState;
};

double TThetaEstimatorData::fisherInfo(const TGenotypeProbabilities &pGenotype, const TGenotypeData deriv_pGenotype) {
	// sum Ri over all sites
	double fInfo = 0.0;
	_begin();
	do {
		// calc Ri
		const double Ri_a = weightedSum(_GL(), deriv_pGenotype);
		const double Ri_b = weightedSum(_GL(), pGenotype);
		const double Ri   = Ri_a / Ri_b;

		// add to Fisher Info
		fInfo += Ri * (Ri + 1.0);
	} while (_next());

	return fInfo;
}

//-------------------------------------------------------
// TThetaEstimatorDataVector
//-------------------------------------------------------
TThetaEstimatorDataVector::TThetaEstimatorDataVector() : TThetaEstimatorData(){};

void TThetaEstimatorDataVector::_saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) {
	// store emission probabilities
	_sites.push_back(genoLik);
};

void TThetaEstimatorDataVector::_emptyStorage() { _sites.clear(); };

void TThetaEstimatorDataVector::_readNext() {
	++_curSite;

	if (_curSite >= _sites.size()) _readState = false;
}

void TThetaEstimatorDataVector::_start() {
	_readState = true;
};

GenotypeLikelihoods::TGenotypeLikelihoods &TThetaEstimatorDataVector::_GL() { return _sites[_curSite]; }

//-------------------------------------------------------
// TThetaEstimatorDataFile
//-------------------------------------------------------
	TThetaEstimatorDataFile::TThetaEstimatorDataFile(std::string TmpFileName) : TThetaEstimatorData(), _sites(TmpFileName) {
	_dataFileName = TmpFileName;
	_sites.openForWriting();
};

void TThetaEstimatorDataFile::_emptyStorage() { _sites.clean(); };

void TThetaEstimatorDataFile::_saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) {
	_sites.save(genoLik);
};

void TThetaEstimatorDataFile::_readNext() {
	++_curSite;
	if (_curSite >= sizeWithData())
		_readState = false;
	else
		_readState = _sites.read(_genotypeLikelihoods);
}

void TThetaEstimatorDataFile::_start() {
	_sites.openForReading();
	--_curSite;
	_readNext(); // read first! This is required to match begin() of a vector
}

GenotypeLikelihoods::TGenotypeLikelihoods &TThetaEstimatorDataFile::_GL() {
	return _genotypeLikelihoods;
};

}; // namespace GenotypeLikelihoods
