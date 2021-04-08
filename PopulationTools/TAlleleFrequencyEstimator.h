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
#include "../TQualityMap.h"
#include "TRandomGenerator.h"

//------------------------------------------------
//THardyWeinbergGenotypeProbabilities
//------------------------------------------------
class THardyWeinbergGenotypeProbabilities{
public:
	double genotypeProbabilities[3];
	double f;
	double oneMinusf;

	THardyWeinbergGenotypeProbabilities(){
		set(0.5);
	};

	void set(double F){
		f = F;
		if(f < 0.0)      f = 0.0;
		else if(f > 1.0) f = 1.0;
		oneMinusf = 1.0 - f;
		genotypeProbabilities[0] = oneMinusf * oneMinusf;
		genotypeProbabilities[1] = 2.0 * f * oneMinusf;
		genotypeProbabilities[2] = f * f;
	};

	double operator[](int g){
		return genotypeProbabilities[g];
	};

	double calcLogLikelihood(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, TGlfConverter & glfConverter);
};

//------------------------------------------------
//TAlleleHardyWeinbergFreqEstimator
//------------------------------------------------
class TAlleleFreqEstimatorHardyWeinberg{
private:
	THardyWeinbergGenotypeProbabilities pGenotype;
	int maxIter;

public:
	TAlleleFreqEstimatorHardyWeinberg();
	double estimate(const TSampleLikelihoods* storage, const uint32_t numSamplesInPop, double epsilonF, TGlfConverter & glfConverter);
};

//------------------------------------------------
//TAlleleFreqEstimatorBayes
//------------------------------------------------
class TAlleleFreqEstimatorBayes{
private:
	TRandomGenerator* randomGenerator;
	TGlfConverter* glfConverter;

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
	THardyWeinbergGenotypeProbabilities pGenotype;
	double f_MAP;
	double logDensity_atMAP;
	double f_CI_lower, f_CI_upper;
	double* f_initialGrid;
	double* f_grid;
	double* density_initialGrid;
	double* density_grid;
	double minPriorSupport, maxPriorSupport;
	double priorDensAtMin, priorDensAtMax;
	double* mcmcSamples;
	bool mcmcSamplesInitialized;

	double _guessInitialAlleleFrequency(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation);
	double _prior(const double & f);
	double _prior(const THardyWeinbergGenotypeProbabilities & pGenotype);
	double _calcPosterior(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, THardyWeinbergGenotypeProbabilities & pGenotype);
	void _fillInitialGrid(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation);
	void _estimateMAP(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation);
	void _estimateCredibleIntervals(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation);

public:
	TAlleleFreqEstimatorBayes(TParameters & Parameters, TLog* logfile, TGlfConverter* GlfConverter, TRandomGenerator* RandomGenerator);
	~TAlleleFreqEstimatorBayes();
	double estimate(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation);
	void composeHeader(std::vector<std::string> & header, const std::string & popName);
	void estimateAndWrite(const TSampleLikelihoods* storage, const uint32_t numSamplesInPop, TOutputFile & out);

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
	TOutputFile outFile;

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
	TLog* logfile;

	// data on individuals
	TPopulationLikelihoodReaderLocus reader;
	TPopulationSamples samples;
	TPopulationLikehoodLocus storage;
	TGlfConverter glfConverter;

	//to measure run time
	struct timeval start;

	//data on loci
	//long _numLoci;

	void _openVCF(TParameters & Parameters);
	void _closeVCF();
	std::vector<std::string> _composeHeaderAlleleFreq(bool writeGenoFreq, bool doBayesian, TAlleleFreqEstimatorBayes* BHWEstimator);
	void _writeBayesianEstimatesOnePop(TOutputFile & out, TSampleLikelihoods* samples, uint32_t numSamples, TAlleleFreqEstimatorBayes* BHWEstimator);
	void _writeEstimatesOnePop(TOutputFile & out, TGenotypeFrequencies & genoFrequencies, double alleleFrequency, TSampleLikelihoods* samples, uint32_t numSamples, TAlleleFreqEstimatorHardyWeinberg & MLHWEstimator, TAlleleFreqEstimatorBayes* BHWEstimator, double epsF, bool writeGenoFreq, bool doBayesian);
	std::vector<std::string> _composeHeaderAlleleFreqComparison(TAlleleFreqEstimatorBayes & BHWEstimator);

public:
	TAlleleFreqEstimator(TParameters & Parameters, TLog* logfile);
	void estimateAlleleFreq(TParameters & Parameters, TRandomGenerator* randomGenerator);
	void compareAlleleFreq(TParameters & Parameters, TRandomGenerator* randomGenerator);
	void writeAlleleFrequencyLikelihoods(TParameters & Parameters, TRandomGenerator* randomGenerator);
};



#endif /* POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_ */
