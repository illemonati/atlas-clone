/*
 * TAlleleFrequencyEstimator.h
 *
 *  Created on: May 13, 2019
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_
#define POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_


#include "TParameters.h"
#include "TPopulationLikelihoods.h"
#include "TRandomGenerator.h"
#include "THardyWeinbergGenotypeProbabilities.h"

namespace PopulationTools{

using coretools::TParameters;
using coretools::TLog;
using coretools::TRandomGenerator;
using coretools::Probability;
using coretools::TOutputFile;

//------------------------------------------------
//TAlleleFreqEstimatorHardyWeinberg
//------------------------------------------------
class TAlleleFreqEstimatorHardyWeinberg{
private:
	uint32_t maxIter;

public:
	TAlleleFreqEstimatorHardyWeinberg();
	Probability estimate(const TSampleLikelihoods* storage, const uint32_t & numSamplesInPop, const double & epsilonF);
};

//------------------------------------------------
//TAlleleFreqEstimatorBayes
//------------------------------------------------

struct FrequencyGridPoint{
	Probability f;
	double density;
};

class TAlleleFreqEstimatorBayes{
private:
	TRandomGenerator* randomGenerator;

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
	Probability f_MAP;
	double logDensity_atMAP;
	double f_CI_lower, f_CI_upper;
	std::vector<FrequencyGridPoint> _initialGrid;
	std::vector<FrequencyGridPoint> _grid;
	double minPriorSupport, maxPriorSupport;
	double priorDensAtMin, priorDensAtMax;

	Probability _guessInitialAlleleFrequency(const TSampleLikelihoods* storage, const uint32_t & numSamplesInPopulation);
	double _prior(const Probability & f) const;
	double _prior(const genometools::THardyWeinbergGenotypeProbabilities & pGenotype) const;
	coretools::LogProbability _calcLogLikelihood(const TSampleLikelihoods* storage, const uint32_t & numSamplesInPopulation, const genometools::THardyWeinbergGenotypeProbabilities & pGenotype);
	double _calcPosterior(const TSampleLikelihoods* storage, const uint32_t & numSamplesInPopulation, const genometools::THardyWeinbergGenotypeProbabilities & pGenotype);
	void _fillInitialGrid(const TSampleLikelihoods* storage, const uint32_t & numSamplesInPopulation);
	void _estimateMAP(const TSampleLikelihoods* storage, const uint32_t & numSamplesInPopulation);
	void _estimateCredibleIntervals(const TSampleLikelihoods* storage, const uint32_t & numSamplesInPopulation);

public:
	TAlleleFreqEstimatorBayes(TParameters & Parameters, TLog* logfile, TRandomGenerator* RandomGenerator);
	~TAlleleFreqEstimatorBayes() = default;
	Probability estimate(const TSampleLikelihoods* storage, const uint32_t & numSamplesInPopulation);
	void composeHeader(std::vector<std::string> & header, const std::string & popName);
	void estimateAndWrite(const TSampleLikelihoods* storage, const uint32_t & numSamplesInPop, TOutputFile & out);

	double credibleIntervalUsed(){ return credibleInterval; };
	double MAP(){ return f_MAP; };
	double lowerCredibleInterval(){ return f_CI_lower; };
	double upperCredibleInterval(){ return f_CI_upper; };
	double runMCMC(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, const double frac, std::vector<double> & mcmcSamples);
	double calcPosteriorf1smallerf2(std::vector<double> & mcmc1, std::vector<double> & mcmc2);
};

//------------------------------------------------
//TAlleleFreqMCMCOutput
//------------------------------------------------
class TAlleleFreqMCMCOutput{
private:
	std::vector<uint32_t> popIndex;
	std::vector<std::string> header;
	std::string outputName;
	coretools::TOutputFile outFile;

public:
	TAlleleFreqMCMCOutput(){}
	TAlleleFreqMCMCOutput(std::string popString, TPopulationSamples & samples, std::string OutputName, TLog* logfile){
		initialize(popString, samples, OutputName, logfile);
	};
	void initialize(std::string popString, TPopulationSamples & samples, std::string OutputName, TLog* logfile);
	void write(std::vector< std::vector<double> > & mcmc, const std::string chr, const uint64_t pos);
};

//------------------------------------------------
//TAlleleFreqEstimator
//------------------------------------------------
class TAlleleFreqEstimator{
private:
	// about vcf-file
	std::string vcfFilename;
	bool vcfRead;
	coretools::TLog* logfile;

	// data on individuals
	TPopulationLikelihoodReaderLocus reader;
	TPopulationSamples samples;
	TPopulationLikehoodLocus storage;

	void _openVCF(TParameters & Parameters);
	void _closeVCF();
	std::vector<std::string> _composeHeaderAlleleFreq(bool writeGenoFreq, bool doBayesian, TAlleleFreqEstimatorBayes* BHWEstimator);
	void _writeBayesianEstimatesOnePop(TOutputFile & out, TSampleLikelihoods* theseSamples, const uint32_t & numSamples, TAlleleFreqEstimatorBayes* BHWEstimator);
	void _writeEstimatesOnePop(TOutputFile & out, genometools::TGenotypeFrequencies & genoFrequencies, double alleleFrequency, TSampleLikelihoods* theseSamples, const uint32_t & numSamples, TAlleleFreqEstimatorHardyWeinberg & MLHWEstimator, TAlleleFreqEstimatorBayes* BHWEstimator, double epsF, bool writeGenoFreq, bool doBayesian);
	std::vector<std::string> _composeHeaderAlleleFreqComparison(TAlleleFreqEstimatorBayes & BHWEstimator);

public:
	TAlleleFreqEstimator(TParameters & Parameters, TLog* Logfile);
	void estimateAlleleFreq(TParameters & Parameters, TRandomGenerator* randomGenerator);
	void compareAlleleFreq(TParameters & Parameters, TRandomGenerator* randomGenerator);
	void writeAlleleFrequencyLikelihoods(TParameters & Parameters, TRandomGenerator* randomGenerator);
};

//--------------------------------------
// Tasks
//--------------------------------------
using coretools::TTask;

class TTask_estimateAlleleFreq:public TTask{
public:
	TTask_estimateAlleleFreq(){ _explanation = "Estimating population allele frequencies"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TAlleleFreqEstimator alleleFreqEstimator(Parameters, Logfile);
		alleleFreqEstimator.estimateAlleleFreq(Parameters, _randomGenerator);
	};
};

class TTask_compareAlleleFreq:public TTask{
public:
	TTask_compareAlleleFreq(){_explanation = "Pairwise comparison of population allele frequencies"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TAlleleFreqEstimator alleleFreqEstimator(Parameters, Logfile);
		alleleFreqEstimator.compareAlleleFreq(Parameters, _randomGenerator);
	};
};

class TTask_alleleFreqLikelihoods:public TTask{
public:
	TTask_alleleFreqLikelihoods(){ _explanation = "Calculating population allele frequency likelihoods under Hardy-Weinberg"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TAlleleFreqEstimator alleleFreqEstimator(Parameters, Logfile);
		alleleFreqEstimator.writeAlleleFrequencyLikelihoods(Parameters, _randomGenerator);
	};
};

}; //end namespace

#endif /* POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_ */
