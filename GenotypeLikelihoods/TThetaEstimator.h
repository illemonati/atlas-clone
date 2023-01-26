/*
 * TThetaEstimator.h
 *
 *  Created on: Sep 25, 2017
 *      Author: phaentu
 */

#ifndef TTHETAESTIMATOR_H_
#define TTHETAESTIMATOR_H_

#include <algorithm>
#include <math.h>
#include <ostream>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <vector>

#include "genometools/GenotypeTypes.h"
#include "coretools/Files/TOutputFile.h"
#include "TGenotypeData.h"
#include "TThetaEstimatorData.h"
#include "TWindow.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenotypeLikelihoods {
class TGenotypeLikelihoodCalculator;
}
namespace GenotypeLikelihoods {
class TSite;
}

namespace GenotypeLikelihoods {

//---------------------------------------------------------------
// Theta
//---------------------------------------------------------------
struct Theta {
	double theta, expMinusTheta, logTheta, thetaConfidence, LL;
	TBaseProbabilities baseFreq;

	Theta() {
		theta           = 0.0;
		thetaConfidence = 0.0;
		expMinusTheta        = 0.0;
		logTheta        = -9e100;
		LL              = -9e100;
	};

	void setTheta(const double val) {
		theta    = val;
		expMinusTheta = exp(-theta);
		logTheta = log(theta);
		LL       = -9e100;
	};

	void setExpTheta(const double val) {
		expMinusTheta = val;
		theta    = -log(val);
		logTheta = log(theta);
		LL       = -9e100;
	};

	void setLogTheta(const double val) {
		logTheta = val;
		theta    = exp(val);
		expMinusTheta = exp(-theta);
		LL       = -9e100;
	};

	void setLogTheta(double val, double newLL) {
		logTheta = val;
		theta    = exp(val);
		expMinusTheta = exp(-theta);
		LL       = newLL;
	};

	std::string getBaseFrequencyString() {
		using genometools::Base;
		return coretools::str::toString("Pi(A) = ", baseFreq[Base::A], ", Pi(C) = ", baseFreq[Base::C],
										", Pi(G) = ", baseFreq[Base::G], ", Pi(T) = ", baseFreq[Base::T]);
	};
};

//---------------------------------------------------------------
// TThetaEstimator
//---------------------------------------------------------------

GenotypeLikelihoods::TGenotypeProbabilities getPGenotype(double expTheta, const TBaseProbabilities &baseFrequencies);
GenotypeLikelihoods::TGenotypeProbabilities getPGenotype(const Theta &thisTheta);

class TThetaEstimator_base {
protected:
	// data
	TThetaEstimatorData *data = nullptr;
	bool dataInitialized      = false;
	bool useTmpFile;
	std::string tmpFileName;

	// initial theta
	double initialTheta              = 0.01;
	double initThetaSearchFactor     = 10;
	int initThetaNumSearchIterations = 100;

	// estimation
	int minSitesWithData;
	GenotypeLikelihoods::TGenotypeProbabilities _pGenotype; // P(g|pi, theta)
	Theta theta;
	bool extraVerbose;

	void initDataStorage();
	void readParametersRegardingInitialSearch();

	void findGoodStartingTheta(TThetaEstimatorData *thisData, Theta &thisTheta, const std::string &tag);

public:
	TThetaEstimator_base();
	TThetaEstimator_base(const TThetaEstimator_base &other);

	virtual ~TThetaEstimator_base() {
		if (dataInitialized) delete data;
	};

	TThetaEstimatorData *pointerToDataContainer() { return data; };

	GenotypeLikelihoods::TGenotypeProbabilities pGenotype() { return ::GenotypeLikelihoods::getPGenotype(theta); };
};

//---------------------------------------------------------------
// TThetaEstimator
//---------------------------------------------------------------
class TThetaEstimator : public TThetaEstimator_base {
private:
	// EM parameters
	int numIterations;
	int numThetaOnlyUpdates;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
	bool estimationSuccessful;
	double _expectedHet;

	// tmp vectors
	GenotypeLikelihoods::TGenotypeData P_G; // see paper

	double _calcFisherInfo(const TGenotypeProbabilities &pGenotype, const TGenotypeData deriv_pGenotype);
	bool _NRAllParams();
	void _NROnlyTheta();
	void _runEMForTheta();
	void _estimateConfidenceInterval();
	void _calcExpectedHet();

public:
	TThetaEstimator();
	TThetaEstimator(const TThetaEstimator &other);

	virtual ~TThetaEstimator(){};

	void clear();
	void add(const TSite &site, const GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods);
	void add(const TWindow &window, const TGenotypeLikelihoodCalculator &glCalculator);
	long sizeWithData() { return data->sizeWithData(); };
	bool estimateTheta();
	void setTheta(const double Theta);
	void setBaseFreq(const GenotypeLikelihoods::TBaseProbabilities &BaseFreq);
	void addToHeader(std::vector<std::string> &header, const std::string &prefix = "");
	void writeEstimateFrequenciesAndTheta(coretools::TOutputFile &out);
	void writeResultsToFile(coretools::TOutputFile &out);
	void calcLikelihoodSurface(coretools::TOutputFile &out, uint32_t steps);
	void bootstrapTheta();
};

//---------------------------------------------------------------
// TThetaEstimatorRatio
//---------------------------------------------------------------
class TThetaEstimatorRatio : public TThetaEstimator_base {
private:
	// second data
	TThetaEstimatorData *data2;
	bool data2Initialized;
	Theta theta2;

	// MCMC parameters
	double phiPriorMean;
	double phiPriorVar;
	double phiPriorOneOverTwoVar;
	double sdProposalKernelTheta1;
	double sdProposalKernelTheta2;
	double sdProposalKernelBaseFreq1;
	double sdProposalKernelBaseFreq2;
	int numIterations;
	int burnin;
	int thinning;
	int numAcceptedTheta1;
	int numAcceptedTheta2;
	int numAcceptedBaseFreq1;
	int numAcceptedBaseFreq2;

	void clearCounters();
	void concludeAcceptanceRate(int numAccepted, int length, const std::string &name);
	void concludeAcceptanceRateUpdateProposal(int numAccepted, int length, double &sd, const std::string &name);
	void concludeAcceptanceRates(int length);
	void concludeAcceptanceRatesUpdateProposal(int length);
	bool updateTheta(TThetaEstimatorData *thisData, Theta &thisTheta, double otherLogThetaMean,
					 double thisSdProposalKernel);
	bool updateBaseFrequencies(TThetaEstimatorData *thisData, Theta &thisTheta, double thisSdProposalKernel);
	void oneMCMCIteration();

public:
	TThetaEstimatorRatio();
	~TThetaEstimatorRatio() {
		if (data2Initialized) delete data2;
	};

	TThetaEstimatorData *pointerToDataContainer2() { return data2; };

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
		std::vector<std::string> header = {"Chr", "start", "end"};

		// add headers of all estimators
		for (size_t i = 0; i < _thetaEstimators.size(); ++i) { _thetaEstimators[i]->addToHeader(header, _prefixes[i]); }
		_out.writeHeader(header);
	};

	void _writeEstimates() {
		for (TThetaEstimator *est : _thetaEstimators) { est->writeResultsToFile(_out); }
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

	void write(const TWindow_base &window) {
		_out.write(window.chrName(), window.from().position() + 1, window.to().position());
		_writeEstimates();
	};

	void write(const std::string &chr, const std::string &start, const std::string &end) {
		_out << chr << start << end;
		_writeEstimates();
	};
};

}; // namespace GenotypeLikelihoods

#endif /* TTHETAESTIMATOR_H_ */
