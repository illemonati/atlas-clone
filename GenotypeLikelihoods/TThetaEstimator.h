/*
 * TThetaEstimator.h
 *
 *  Created on: Sep 25, 2017
 *      Author: phaentu
 */

#ifndef TTHETAESTIMATOR_H_
#define TTHETAESTIMATOR_H_

#include <memory>
#include <string>
#include <vector>

#include "genometools/Genotypes/Base.h"
#include "coretools/Files/TOutputFile.h"
#include "genometools/Genotypes/Containers.h"
#include "TThetaEstimatorData.h"
#include "TWindow.h"

namespace GenotypeLikelihoods {
class TErrorModels;
}
namespace GenotypeLikelihoods {
class TSite;
}

namespace GenotypeLikelihoods {

//---------------------------------------------------------------
// Theta
//---------------------------------------------------------------
struct Theta {
	double theta           = 0.0;
	double expMinusTheta   = 0.0;
	double logTheta        = -9e100;
	double thetaConfidence = 0.0;
	double LL              = -9e100;
	genometools::TBaseProbabilities baseFreq;

	void set(double val) {
		theta    = val;
		expMinusTheta = exp(-theta);
		logTheta = log(theta);
		LL       = -9e100;
	};

	void setExp(double val) {
		expMinusTheta = val;
		theta    = -log(val);
		logTheta = log(theta);
		LL       = -9e100;
	};

	void setLog(double val) {
		logTheta = val;
		theta    = exp(val);
		expMinusTheta = exp(-theta);
		LL       = -9e100;
	};

	void setLog(double val, double newLL) {
		logTheta = val;
		theta    = exp(val);
		expMinusTheta = exp(-theta);
		LL       = newLL;
	};

	std::string string() {
		using genometools::Base;
		return coretools::str::toString("Pi(A) = ", baseFreq[Base::A], ", Pi(C) = ", baseFreq[Base::C],
										", Pi(G) = ", baseFreq[Base::G], ", Pi(T) = ", baseFreq[Base::T]);
	};
};

//---------------------------------------------------------------
// TThetaEstimator
//---------------------------------------------------------------

genometools::TGenotypeProbabilities getPGenotype(double expTheta, const genometools::TBaseProbabilities &baseFrequencies);
genometools::TGenotypeProbabilities getPGenotype(const Theta &thisTheta);

class TThetaEstimator_base {
protected:
	// data
	std::unique_ptr<TThetaEstimatorData> _data;
	std::string _tmpFileName;
	bool _useTmpFile    = false;
	bool _equalBaseFreqs = false;

	// initial theta
	double _initialTheta              = 0.01;
	double _initThetaSearchFactor     = 10;
	int _initThetaNumSearchIterations = 100;

	// estimation
	size_t _minSitesWithData;
	Theta _theta;
	bool _extraVerbose;

	void _initDataStorage();
	void _readParametersRegardingInitialSearch();
	void _findGoodStartingTheta(TThetaEstimatorData *thisData, Theta &thisTheta, const std::string &tag);

public:
	TThetaEstimator_base();
	TThetaEstimator_base(const TThetaEstimator_base &other);
	virtual ~TThetaEstimator_base() = default;

	TThetaEstimatorData *pointerToDataContainer() { return _data.get(); };
	genometools::TGenotypeProbabilities pGenotype() { return ::GenotypeLikelihoods::getPGenotype(_theta); };
};

//---------------------------------------------------------------
// TThetaEstimator
//---------------------------------------------------------------
class TThetaEstimator : public TThetaEstimator_base {
private:
	// EM parameters
	int _numIterations;
	int _numThetaOnlyUpdates;
	double _maxEpsilon;
	int _NewtonRaphsonNumIterations;
	double _NewtonRaphsonMaxF;
	bool _estimationSuccessful = false;
	double _expectedHet        = 0.0;

	bool _NRAllParams(const genometools::TGenotypeProbabilities &pGenotype);
	void _NROnlyTheta(const genometools::TGenotypeProbabilities &pGenotype);
	void _runEMForTheta();
	void _estimateConfidenceInterval();
	void _calcExpectedHet();

public:
	TThetaEstimator();

	void clear();
	void add(const TSite &site, const genometools::TGenotypeLikelihoods &genotypeLikelihoods);
	void add(const TWindow &window, const TErrorModels &glCalculator);
	long sizeWithData() { return _data->sizeWithData(); };
	bool estimateTheta();
	void setTheta(double Theta);
	void setBaseFreq(const genometools::TBaseProbabilities &BaseFreq);
	void addToHeader(std::vector<std::string> &header, const std::string &prefix = "");
	void writeEstimateFrequenciesAndTheta(coretools::TOutputFile &out);
	void writeResultsToFile(coretools::TOutputFile &out, size_t NumMaskedSites);
	void calcLikelihoodSurface(coretools::TOutputFile &out, size_t steps);
	void bootstrapTheta();
};

//---------------------------------------------------------------
// TThetaEstimatorRatio
//---------------------------------------------------------------
class TThetaEstimatorRatio : public TThetaEstimator_base {
private:
	// second data
	std::unique_ptr<TThetaEstimatorData> _data2;
	Theta _theta2;

	// MCMC parameters
	double _phiPriorMean;
	double _phiPriorVar;
	double _phiPriorOneOverTwoVar;
	double _sdProposalKernelTheta1;
	double _sdProposalKernelTheta2;
	double _sdProposalKernelBaseFreq1;
	double _sdProposalKernelBaseFreq2;
	int _numIterations;
	int _burnin;
	int _thinning;
	int _numAcceptedTheta1;
	int _numAcceptedTheta2;
	int _numAcceptedBaseFreq1;
	int _numAcceptedBaseFreq2;

	void _clearCounters();
	void _concludeAcceptanceRate(int numAccepted, int length, const std::string &name);
	void _concludeAcceptanceRateUpdateProposal(int numAccepted, int length, double &sd, const std::string &name);
	void _concludeAcceptanceRates(int length);
	void _concludeAcceptanceRatesUpdateProposal(int length);
	bool _updateTheta(TThetaEstimatorData *thisData, Theta &thisTheta, double otherLogThetaMean,
					 double thisSdProposalKernel);
	bool _updateBaseFrequencies(TThetaEstimatorData *thisData, Theta &thisTheta, double thisSdProposalKernel);
	void _oneMCMCIteration();

public:
	TThetaEstimatorRatio();

	TThetaEstimatorData *pointerToDataContainer2() { return _data2.get(); };
	void estimateRatio(std::string ouputName);
};

//---------------------------------------------------------------
// TThetaOutputFile
//---------------------------------------------------------------
class TThetaOutputFile {
protected:
	coretools::TOutputFile _out;
	std::vector<TThetaEstimator *> _thetaEstimators;
	std::vector<std::string> _prefixes;

	void _writeHeader() {
		std::vector<std::string> header = {"Chr", "Start", "End"};

		// add headers of all estimators
		for (size_t i = 0; i < _thetaEstimators.size(); ++i) { _thetaEstimators[i]->addToHeader(header, _prefixes[i]); }
		_out.writeHeader(header);
	};

	void _writeEstimates(size_t NumMaskedSites) {
		for (TThetaEstimator *est : _thetaEstimators) { est->writeResultsToFile(_out, NumMaskedSites); }
		_out.endln();
	};

public:
	TThetaOutputFile(){};

	TThetaOutputFile(TThetaEstimator *Estimator, const std::string &Filename) {
		addEstimator(Estimator, "");
		open(Filename);
	};

	~TThetaOutputFile() { close(); };

	coretools::TOutputFile &file() { return _out; };

	void addEstimator(TThetaEstimator *Estimator, const std::string &Prefix) {
		if (_out.isOpen()) { UERROR("Can not add estimators to an open TThetaOutputFile!"); }
		_thetaEstimators.push_back(Estimator);
		_prefixes.push_back(Prefix);
	};

	void open(const std::string &Filename) {
		coretools::instances::logfile().list("Will write theta estimates to file '" + Filename + "'.");
		_out.open(Filename);
		_writeHeader();
		_out.precision(9);
	};

	void open(TThetaEstimator *Estimator, const std::string &Filename) {
		addEstimator(Estimator, "");
		open(Filename);
	};

	void close() { _out.close(); };

	void write(const TWindow &window) {
		_out.write(window.chrName(), window.from().position() + 1, window.to().position());
		_writeEstimates(window.numMaskedSites());
	};

	void write(const std::string &chr, const std::string &start, const std::string &end, size_t TotMaskedSites) {
		_out << chr << start << end;
		_writeEstimates(TotMaskedSites);
	};
};

}; // namespace GenotypeLikelihoods

#endif /* TTHETAESTIMATOR_H_ */
