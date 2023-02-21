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
	gzFile fp = nullptr;
	static constexpr size_t sizeOfData = sizeof(double) * 10;

	bool isOpenForWriting = false;
	bool isOpenForReading = false;
	bool wasWritten = false;

public:
	TThetaEstimatorTemporaryFile(std::string Filename): filename(Filename) {}
	~TThetaEstimatorTemporaryFile() { clean(); };

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
private:
	// counters
	size_t numSitesCoveredTwiceOrMore;
	size_t totNumSitesAdded;
	double cumulativeDepth;
	size_t numSitesWithData;

	std::vector<GenotypeLikelihoods::TBaseData> baseFreqs;
	std::vector<size_t> numBootstrapRepsPerEntry;
	size_t curRep;

protected:
	bool readState;
	size_t curSite;

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
	void clearBootstrap() noexcept {numBootstrapRepsPerEntry.clear();}

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

	void saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) override;
	void readNext() override;
	void emptyStorage() override;
	void _begin() override;

public:
	TThetaEstimatorDataVector();

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
	void _begin() override;

public:
	TThetaEstimatorDataFile(std::string TmpFileName);
	~TThetaEstimatorDataFile() { clear(); };

	GenotypeLikelihoods::TGenotypeLikelihoods &curGenotypeLikelihoods() override;
};

}; // namespace GenotypeLikelihoods

#endif /* TTHETAESTIMATORDATA_H_ */
