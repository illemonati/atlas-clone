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
#include "TVCFReader.h"
#include "TTask.h"

namespace PopulationTools{

//------------------------------------------------
//TAlleleFreqEstimatorHardyWeinberg
//------------------------------------------------
class TAlleleFreqEstimatorHardyWeinberg{
private:
	uint32_t maxIter;

public:
	TAlleleFreqEstimatorHardyWeinberg();
	coretools::Probability estimate(const TSampleLikelihoods* storage, uint32_t numSamplesInPop, double epsilonF);
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
	coretools::TRandomGenerator* randomGenerator;

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

	coretools::Probability _guessInitialAlleleFrequency(const TSampleLikelihoods* storage, uint32_t numSamplesInPopulation);
	double _prior(const coretools::Probability & f) const;
	double _prior(const genometools::THardyWeinbergGenotypeProbabilities & pGenotype) const;
	coretools::LogProbability _calcLogLikelihood(const TSampleLikelihoods* storage, uint32_t numSamplesInPopulation, const genometools::THardyWeinbergGenotypeProbabilities & pGenotype);
	double _calcPosterior(const TSampleLikelihoods* storage, uint32_t numSamplesInPopulation, const genometools::THardyWeinbergGenotypeProbabilities & pGenotype);
	void _fillInitialGrid(const TSampleLikelihoods* storage, uint32_t numSamplesInPopulation);
	void _estimateMAP(const TSampleLikelihoods* storage, uint32_t numSamplesInPopulation);
	void _estimateCredibleIntervals(const TSampleLikelihoods* storage, uint32_t numSamplesInPopulation);

public:
	TAlleleFreqEstimatorBayes(coretools::TParameters & Parameters, coretools::TLog* logfile, coretools::TRandomGenerator* RandomGenerator);
	~TAlleleFreqEstimatorBayes() = default;
	coretools::Probability estimate(const TSampleLikelihoods* storage, uint32_t numSamplesInPopulation);
	void composeHeader(std::vector<std::string> & header, const std::string & popName);
	void estimateAndWrite(const TSampleLikelihoods* storage, uint32_t numSamplesInPop, coretools::TOutputFile & out);

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
	TAlleleFreqMCMCOutput(std::string popString, genometools::TPopulationSamples & samples, std::string OutputName, coretools::TLog* logfile){
		initialize(popString, samples, OutputName, logfile);
	};
	void initialize(std::string popString, genometools::TPopulationSamples & samples, std::string OutputName, coretools::TLog* logfile);
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
    genometools::TPopulationLikelihoodReaderLocus reader;
    genometools::TPopulationSamples samples;
    genometools::TPopulationLikehoodLocus<genometools::HighPrecisionPhredIntProbability> storage;

	void _openVCF(coretools::TParameters & Parameters);
	void _closeVCF();
	std::vector<std::string> _composeHeaderAlleleFreq(bool writeGenoFreq, bool doBayesian, TAlleleFreqEstimatorBayes* BHWEstimator);
	void _writeBayesianEstimatesOnePop(coretools::TOutputFile & out, TSampleLikelihoods* theseSamples, uint32_t numSamples, TAlleleFreqEstimatorBayes* BHWEstimator);
	void _writeEstimatesOnePop(coretools::TOutputFile & out, genometools::TGenotypeFrequencies & genoFrequencies, double alleleFrequency, TSampleLikelihoods* theseSamples, uint32_t numSamples, TAlleleFreqEstimatorHardyWeinberg & MLHWEstimator, TAlleleFreqEstimatorBayes* BHWEstimator, double epsF, bool writeGenoFreq, bool doBayesian);
	std::vector<std::string> _composeHeaderAlleleFreqComparison(TAlleleFreqEstimatorBayes & BHWEstimator);

public:
	TAlleleFreqEstimator(coretools::TParameters & Parameters, coretools::TLog* Logfile);
	void estimateAlleleFreq(coretools::TParameters & Parameters, coretools::TRandomGenerator* randomGenerator);
	void compareAlleleFreq(coretools::TParameters & Parameters, coretools::TRandomGenerator* randomGenerator);
	void writeAlleleFrequencyLikelihoods(coretools::TParameters & Parameters, coretools::TRandomGenerator* randomGenerator);
};

//--------------------------------------
// Tasks
//--------------------------------------

class TTask_estimateAlleleFreq:public coretools::TTask{
public:
	TTask_estimateAlleleFreq(){ _explanation = "Estimating population allele frequencies"; };

	void run(){
		using namespace coretools::instances;
		TAlleleFreqEstimator alleleFreqEstimator(parameters(), &logfile());
		alleleFreqEstimator.estimateAlleleFreq(parameters(), &randomGenerator());
	};
};

class TTask_compareAlleleFreq:public coretools::TTask{
public:
	TTask_compareAlleleFreq(){_explanation = "Pairwise comparison of population allele frequencies"; };

	void run(){
		using namespace coretools::instances;
		using namespace coretools::instances;
		TAlleleFreqEstimator alleleFreqEstimator(parameters(), &logfile());
		alleleFreqEstimator.compareAlleleFreq(parameters(), &randomGenerator());
	};
};

class TTask_alleleFreqLikelihoods:public coretools::TTask{
public:
	TTask_alleleFreqLikelihoods(){ _explanation = "Calculating population allele frequency likelihoods under Hardy-Weinberg"; };

	void run(){
		using namespace coretools::instances;
		TAlleleFreqEstimator alleleFreqEstimator(parameters(), &logfile());
		alleleFreqEstimator.writeAlleleFrequencyLikelihoods(parameters(), &randomGenerator());
	};
};

}; //end namespace

#endif /* POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_ */
