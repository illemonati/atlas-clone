/*
 * TThetaEstimatorData.h
 *
 *  Created on: Jun 17, 2018
 *      Author: phaentu
 */

#ifndef TTHETAESTIMATORDATA_H_
#define TTHETAESTIMATORDATA_H_

#include <string>
#include <vector>
#include <zlib.h>

#include "GenotypeData.h"
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
	static constexpr size_t _sizeOfData = sizeof(double) * 10;

	std::string _filename;
	gzFile _fp             = nullptr;
	bool _isOpenForWriting = false;
	bool _isOpenForReading = false;
	bool _wasWritten       = false;

public:
	TThetaEstimatorTemporaryFile(std::string Filename) : _filename(Filename) {}
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
	size_t _numSites;
	size_t _numSitesData;
	size_t _numSites2x;
	double _cumulativeDepth;

	std::vector<GenotypeLikelihoods::TBaseData> _baseFreqs;
	std::vector<size_t> _numBootstrapRepsPerEntry;
	size_t _curRep;

protected:
	bool _readState;
	size_t _curSite;

	bool _isBootstrapped() const noexcept { return !_numBootstrapRepsPerEntry.empty(); }
	bool _begin();
	bool _next();

	virtual void _saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &) = 0;
	virtual void _start()                                                     = 0;
	virtual void _readNext()                                                  = 0;
	virtual void _emptyStorage()                                              = 0;
	virtual GenotypeLikelihoods::TGenotypeLikelihoods &_GL()                  = 0;

public:
	TThetaEstimatorData();
	virtual ~TThetaEstimatorData() = default;

	void add(const GenotypeLikelihoods::TSite &site, const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik);
	void clear();

	void bootstrap();
	void clearBootstrap() noexcept { _numBootstrapRepsPerEntry.clear(); }

	size_t size() { return _numSites; };
	size_t sizeWithData() { return _numSitesData; };

	void addToHeader(std::vector<std::string> &header, const std::string &prefix);
	void writeSite(coretools::TOutputFile &out, size_t NumMaskedSites);
	TBaseProbabilities baseFrequencies();

	double fisherInfo(const TGenotypeProbabilities &pGenotype, const TGenotypeData deriv_pGenotype);


	GenotypeLikelihoods::TGenotypeData P_G(const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype);
	double calcLogLikelihood(const GenotypeLikelihoods::TGenotypeProbabilities &pGenotype);
};

//---------------------------------------------------------------
// TThetaEstimatorDataVector
//---------------------------------------------------------------
class TThetaEstimatorDataVector final : public TThetaEstimatorData {
private:
	std::vector<GenotypeLikelihoods::TGenotypeLikelihoods> _sites;

	void _saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) override;
	void _readNext() override;
	void _emptyStorage() override;
	void _start() override;
	GenotypeLikelihoods::TGenotypeLikelihoods &_GL() override;

public:
	TThetaEstimatorDataVector();

};

//---------------------------------------------------------------
// TThetaEstimatorDataFile
//---------------------------------------------------------------
class TThetaEstimatorDataFile final : public TThetaEstimatorData {
protected:
	std::string _dataFileName;

	TThetaEstimatorTemporaryFile _sites;
	GenotypeLikelihoods::TGenotypeLikelihoods _genotypeLikelihoods;

	void _saveSite(const GenotypeLikelihoods::TGenotypeLikelihoods &genoLik) override;
	void _readNext() override;
	void _emptyStorage() override;
	void _start() override;
	GenotypeLikelihoods::TGenotypeLikelihoods &_GL() override;

public:
	TThetaEstimatorDataFile(std::string TmpFileName);
	~TThetaEstimatorDataFile() { clear(); };

};

}; // namespace GenotypeLikelihoods

#endif /* TTHETAESTIMATORDATA_H_ */
