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
	long numSitesCoveredTwiceOrMore;
	long totNumSitesAdded;
	double cumulativeDepth;

	GenotypeLikelihoods::TBaseData tmpBaseFreq;
	bool isBootstrapped;
	long numBootstrappedSites;
	int maxKforPoissonPlusOne;
	double *poissonProb;
	uint8_t *numBootstrapRepsPerEntry;
	bool numBootstrapRepsPerEntryInitialized;

	bool readState;
	long curSite;
	uint8_t curRep;

	// tmp variables
	int g;
	double sum;

	virtual void saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &) {
		throw "Not available in TThetaEstimatorData base class!";
	};
	virtual void emptyStorage(){};
	void fillPoissonForBootstrap(const double lambda);
	virtual void _begin() { throw "Not available in TThetaEstimatorData base class!"; };
	virtual void readNext() { throw "Not available in TThetaEstimatorData base class!"; };

public:
	TThetaEstimatorData();
	virtual ~TThetaEstimatorData() {
		clear();
		delete[] poissonProb;
		if (numBootstrapRepsPerEntryInitialized) delete[] numBootstrapRepsPerEntry;
	};
	long numSitesWithData;

	void add(const GenotypeLikelihoods::TSite &site, const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik);
	void clear();

	void bootstrap();
	void clearBootstrap();

	virtual bool begin();
	virtual bool next();
	virtual bool isEnd() { throw "Not available in TThetaEstimatorData base class!"; };
	virtual GenotypeLikelihoods::TGenotypeLikelihoods &curGenotypeLikelihoods() {
		throw "Not available in TThetaEstimatorData base class!";
	};

	long size() { return totNumSitesAdded; };
	long sizeWithData() {
		if (isBootstrapped) return numBootstrappedSites;
		return numSitesWithData;
	};

	void addToHeader(std::vector<std::string> &header, const std::string prefix);
	void writeSite(coretools::TOutputFile &out);
	TBaseProbabilities baseFrequencies();
	virtual void fillP_G(GenotypeLikelihoods::TGenotypeData &P_G,
						 const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype);
	virtual double calcLogLikelihood(const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype);
};

//---------------------------------------------------------------
// TThetaEstimatorDataVector
//---------------------------------------------------------------
class TThetaEstimatorDataVector : public TThetaEstimatorData {
private:
	std::vector<GenotypeLikelihoods::TGenotypeLikelihoods> sites;
	std::vector<GenotypeLikelihoods::TGenotypeLikelihoods>::iterator siteIt;

	void saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) override;
	void emptyStorage() override;
	void readNext() override;

public:
	TThetaEstimatorDataVector();
	virtual ~TThetaEstimatorDataVector() { clear(); };

	void _begin() override;
	bool isEnd() override;
	GenotypeLikelihoods::TGenotypeLikelihoods &curGenotypeLikelihoods() override;
	void fillP_G(GenotypeLikelihoods::TGenotypeData &P_G,
				 const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype) override;
	double calcLogLikelihood(const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype) override;
};

//---------------------------------------------------------------
// TThetaEstimatorDataFile
//---------------------------------------------------------------
class TThetaEstimatorDataFile : public TThetaEstimatorData {
protected:
	std::string dataFileName;

	TThetaEstimatorTemporaryFile sites;
	GenotypeLikelihoods::TGenotypeLikelihoods genotypeLikelihoods;

	void saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) override;
	void emptyStorage() override;
	void readNext() override;

public:
	TThetaEstimatorDataFile(std::string TmpFileName);
	~TThetaEstimatorDataFile() { clear(); };

	void _begin() override;
	bool isEnd() override;
	GenotypeLikelihoods::TGenotypeLikelihoods &curGenotypeLikelihoods() override;
};

}; // namespace GenotypeLikelihoods

#endif /* TTHETAESTIMATORDATA_H_ */
