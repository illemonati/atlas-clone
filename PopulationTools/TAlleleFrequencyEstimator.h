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
#include "TQualityMap.h"

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

	double guessInitialAlleleFrequency(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, TGlfConverter & glfConverter);
	double _prior(const double & f);
	double _prior(const THardyWeinbergGenotypeProbabilities & pGenotype);
	double calcPosterior(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, THardyWeinbergGenotypeProbabilities & pGenotype, TGlfConverter & glfConverter);
	void fillInitialGrid(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, TGlfConverter & glfConverter);
	void estimateMAP(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, TGlfConverter & glfConverter);
	void estimateCredibleIntervals(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, TGlfConverter & glfConverter);

public:
	TAlleleFreqEstimatorBayes(TParameters & Parameters, TLog* logfile, TRandomGenerator* RandomGenerator);
	~TAlleleFreqEstimatorBayes();
	double estimate(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, TGlfConverter & glfConverter);

	double credibleIntervalUsed(){ return credibleInterval; };
	double MAP(){ return f_MAP; };
	double lowerCredibleInterval(){ return f_CI_lower; };
	double upperCredibleInterval(){ return f_CI_upper; };
	double runMCMC(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, TGlfConverter & glfConverter, const uint32_t numIterations, const double frac);
	double* getMcmcSamples(){ return mcmcSamples; };
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
	TPopulationSamples samples;
	TGlfConverter glfConverter;

	//data on loci
	long _numLoci;

	std::vector<std::string> writeHeaderAlleleFreq(bool writeGenoFreq, bool doBayesian, TAlleleFreqEstimatorBayes* BHWEstimator);
	void writeEstimatesOnePop(TOutputFileZipped & out, TGenotypeFrequencies & genoFrequencies, double alleleFrequency, TSampleLikelihoods* samples, int numSamples, TAlleleFreqEstimatorHardyWeinberg & MLHWEstimator, TAlleleFreqEstimatorBayes* BHWEstimator, double epsF, bool writeGenoFreq, bool doBayesian);
	std::vector<std::string> writeHeaderAlleleFreqComparison(bool writeGenoFreq, TAlleleFreqEstimatorBayes* BHWEstimator);

public:
	TAlleleFreqEstimator(TParameters & Parameters, TLog* logfile);
	void estimateAlleleFreq(TParameters & Parameters, TRandomGenerator* randomGenerator);
	void compareAlleleFreq(TParameters & Parameters, TRandomGenerator* randomGenerator);
};



#endif /* POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_ */
