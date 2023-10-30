/*
 * TAlleleFrequencyEstimator.h
 *
 *  Created on: May 13, 2019
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_
#define POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_

#include <stdint.h>
#include <string>
#include <vector>

#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TTask.h"
#include "coretools/Types/probability.h"
#include "genometools/PhredProbabilityTypes.h"
#include "genometools/THardyWeinbergGenotypeProbabilities.h"
#include "genometools/TSampleLikelihoods.h"
#include "genometools/VCF/TPopulation.h"
#include "genometools/VCF/TPopulationLikelihoodLocus.h"
#include "genometools/VCF/TPopulationLikelihoods.h"

namespace genometools { class TGenotypeFrequencies; }


using TSampleLikelihoods = genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability>;

namespace PopulationTools{

//------------------------------------------------
//TAlleleFreqEstimatorHardyWeinberg
//------------------------------------------------
class TAlleleFreqEstimatorHardyWeinberg{
private:
	size_t maxIter;
	coretools::Probability _alleleFrequency;

public:
	TAlleleFreqEstimatorHardyWeinberg();
	coretools::Probability estimate(const TSampleLikelihoods* storage, size_t numSamplesInPop, double epsilonF);
	coretools::Probability alleleFrequency(){ return _alleleFrequency; }
	coretools::Log10Probability calculateLog10Likelihood(const TSampleLikelihoods* storage, size_t numSamplesInPop) const noexcept;
};

//------------------------------------------------
//TAlleleFreqEstimatorBayes
//------------------------------------------------

struct FrequencyGridPoint{
	coretools::Probability f;
	double density;
};

class TAlleleFreqEstimatorBayes{
private:

	double alpha, beta;
	double alphaMinusOne, betaMinusOne;

	//MAP and CI search
	int numMAPSIterations;
	int initialGridSize;
	int initialGridLast;
	double logGridThreshold;
	int gridSize;
	int gridLast;
	double credibleInterval;
	coretools::Probability f_MAP;
	double logDensity_atMAP;
	double f_CI_lower, f_CI_upper;
	std::vector<FrequencyGridPoint> _initialGrid;
	std::vector<FrequencyGridPoint> _grid;
	double minPriorSupport, maxPriorSupport;
	double priorDensAtMin, priorDensAtMax;

	coretools::Probability _guessInitialAlleleFrequency(const TSampleLikelihoods* storage, size_t numSamplesInPopulation);
	double _prior(const coretools::Probability & f) const;
	double _prior(const genometools::THardyWeinbergGenotypeProbabilities & pGenotype) const;
	coretools::LogProbability _calcLogLikelihood(const TSampleLikelihoods* storage, size_t numSamplesInPopulation, const genometools::THardyWeinbergGenotypeProbabilities & pGenotype);
	double _calcPosterior(const TSampleLikelihoods* storage, size_t numSamplesInPopulation, const genometools::THardyWeinbergGenotypeProbabilities & pGenotype);
	void _fillInitialGrid(const TSampleLikelihoods* storage, size_t numSamplesInPopulation);
	void _estimateMAP(const TSampleLikelihoods* storage, size_t numSamplesInPopulation);
	void _estimateCredibleIntervals(const TSampleLikelihoods* storage, size_t numSamplesInPopulation);

public:
	TAlleleFreqEstimatorBayes();
	~TAlleleFreqEstimatorBayes() = default;
	coretools::Probability estimate(const TSampleLikelihoods* storage, size_t numSamplesInPopulation);
	void composeHeader(std::vector<std::string> & header, const std::string & popName);
	void estimateAndWrite(const TSampleLikelihoods* storage, size_t numSamplesInPop, coretools::TOutputFile & out);

	double credibleIntervalUsed(){ return credibleInterval; };
	double MAP(){ return f_MAP; };
	double lowerCredibleInterval(){ return f_CI_lower; };
	double upperCredibleInterval(){ return f_CI_upper; };
	double runMCMC(const TSampleLikelihoods* storage, size_t numSamplesInPopulation, double frac, std::vector<double> & mcmcSamples);
	double calcPosteriorf1smallerf2(std::vector<double> & mcmc1, std::vector<double> & mcmc2);
};

//------------------------------------------------
//TAlleleFreqMCMCOutput
//------------------------------------------------
class TAlleleFreqMCMCOutput{
private:
	std::vector<size_t> popIndex;
	std::vector<std::string> header;
	std::string outputName;
	coretools::TOutputFile outFile;

public:
	TAlleleFreqMCMCOutput(){}
	TAlleleFreqMCMCOutput(std::string popString, genometools::TPopulationSamples & samples, std::string OutputName){
		initialize(popString, samples, OutputName);
	};
	void initialize(std::string popString, genometools::TPopulationSamples & samples, std::string OutputName);
	void write(std::vector< std::vector<double> > & mcmc, std::string_view chr, size_t pos);
};

//------------------------------------------------
//TAlleleFreqEstimator
//------------------------------------------------
class TAlleleFreqEstimator{
private:
	// about vcf-file
	std::string vcfFilename;
	bool vcfRead;

	// data on individuals
    genometools::TPopulationLikelihoodReaderLocus reader;
    genometools::TPopulationSamples samples;
    genometools::TPopulationLikehoodLocus<TSampleLikelihoods> storage;

	void _openVCF();
	void _closeVCF();
	std::vector<std::string> _composeHeaderAlleleFreq(bool writeGenoFreq, bool doBayesian, TAlleleFreqEstimatorBayes* BHWEstimator, bool writeLikelihoods);
	void _writeBayesianEstimatesOnePop(coretools::TOutputFile & out, TSampleLikelihoods* theseSamples, size_t numSamples, TAlleleFreqEstimatorBayes* BHWEstimator);
	void _writeEstimatesOnePop(coretools::TOutputFile & out, genometools::TGenotypeFrequencies & genoFrequencies, TSampleLikelihoods* theseSamples, size_t numSamples, TAlleleFreqEstimatorHardyWeinberg & MLHWEstimator, TAlleleFreqEstimatorBayes* BHWEstimator, double epsF, bool writeGenoFreq, bool doBayesian, bool writeLikelihoods);
	std::vector<std::string> _composeHeaderAlleleFreqComparison(TAlleleFreqEstimatorBayes & BHWEstimator);

public:
	TAlleleFreqEstimator();
	void estimateAlleleFreq();
	void compareAlleleFreq();
	void writeAlleleFrequencyLikelihoods();
	void run() {
		 if (coretools::instances::parameters().exists("compare")) {
			compareAlleleFreq();
		} else if  (coretools::instances::parameters().exists("likelihoods")) {
			writeAlleleFrequencyLikelihoods();
		} else {
			estimateAlleleFreq();
		}
	}
};

}; //end namespace

#endif /* POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_ */
