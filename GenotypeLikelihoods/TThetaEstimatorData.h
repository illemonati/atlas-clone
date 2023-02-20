/*
 * TThetaEstimatorData.h
 *
 *  Created on: Jun 17, 2018
 *      Author: phaentu
 */

#ifndef TTHETAESTIMATORDATA_H_
#define TTHETAESTIMATORDATA_H_

#include <stdint.h>
#include <string>
#include <vector>
#include <zlib.h>

#include "TGenotypeData.h"
#include "coretools/Files/gzstream.h"

namespace GenotypeLikelihoods {
class TSite;
}
namespace coretools {
class TOutputFile;
}

namespace GenotypeLikelihoods {

//---------------------------------------------------------------
// TThetaEstimatorTemporaryFile
//---------------------------------------------------------------
class TThetaEstimatorTemporaryFile {
private:
	std::string filename;
	gzFile fp;
	int sizeOfData;

	bool isOpen;
	bool isOpenForWriting;
	bool isOpenForReading;
	bool wasWritten;

public:
	TThetaEstimatorTemporaryFile();
	TThetaEstimatorTemporaryFile(std::string filename);
	~TThetaEstimatorTemporaryFile() { clean(); };

	void init(std::string filename);

	void openForWriting();
	void openForReading();
	void close();
	void clean();
	bool isEOF();

	void save(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik);
	bool read(GenotypeLikelihoods::TGenotypeLikelihoods &genoLik);
};

//---------------------------------------------------------------
// TThetaEstimatorData
//---------------------------------------------------------------
class TThetaEstimatorData {
protected:
	// counters
	size_t numSitesCoveredTwiceOrMore;
	size_t totNumSitesAdded;
	double cumulativeDepth;
	size_t numSitesWithData;

	GenotypeLikelihoods::TBaseData tmpBaseFreq;
	std::vector<size_t> numBootstrapRepsPerEntry;

	bool readState;
	size_t curSite;
	size_t curRep;

	bool _isBootstrapped() const noexcept {return !numBootstrapRepsPerEntry.empty();}

	virtual void saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &) = 0;
	virtual void _begin()                                                    = 0;
	virtual void readNext()                                                  = 0;
	virtual void emptyStorage()                                              = 0;

public:
	TThetaEstimatorData();
	virtual ~TThetaEstimatorData() = default;

	void add(const GenotypeLikelihoods::TSite &site, const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik);
	void clear();

	void bootstrap();
	void clearBootstrap();

	bool begin();
	bool next();
	size_t size() { return totNumSitesAdded; };
	size_t sizeWithData() { return numSitesWithData; };

	void addToHeader(std::vector<std::string> &header, const std::string prefix);
	void writeSite(coretools::TOutputFile &out);
	TBaseProbabilities baseFrequencies();

	virtual GenotypeLikelihoods::TGenotypeLikelihoods &curGenotypeLikelihoods() = 0;

	GenotypeLikelihoods::TGenotypeData P_G(const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype);
	double calcLogLikelihood(const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype);
};

//---------------------------------------------------------------
// TThetaEstimatorDataVector
//---------------------------------------------------------------
class TThetaEstimatorDataVector final : public TThetaEstimatorData {
private:
	std::vector<GenotypeLikelihoods::TGenotypeLikelihoods> sites;
	std::vector<GenotypeLikelihoods::TGenotypeLikelihoods>::iterator siteIt;

	void saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) override;
	void readNext() override;
	void emptyStorage() override;

public:
	TThetaEstimatorDataVector();

	void _begin() override;
	GenotypeLikelihoods::TGenotypeLikelihoods &curGenotypeLikelihoods() override;
};

//---------------------------------------------------------------
// TThetaEstimatorDataFile
//---------------------------------------------------------------
class TThetaEstimatorDataFile final : public TThetaEstimatorData {
protected:
	std::string dataFileName;

	TThetaEstimatorTemporaryFile sites;
	GenotypeLikelihoods::TGenotypeLikelihoods genotypeLikelihoods;

	void saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) override;
	void readNext() override;
	void emptyStorage() override;

public:
	TThetaEstimatorDataFile(std::string TmpFileName);
	~TThetaEstimatorDataFile() { clear(); };

	void _begin() override;
	GenotypeLikelihoods::TGenotypeLikelihoods &curGenotypeLikelihoods() override;
};

}; // namespace GenotypeLikelihoods

#endif /* TTHETAESTIMATORDATA_H_ */
