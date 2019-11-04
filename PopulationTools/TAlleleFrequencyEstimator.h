/*
 * TAlleleFrequencyEstimator.h
 *
 *  Created on: May 13, 2019
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_
#define POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_


#include "../TParameters.h"
#include "TPopulationLikelihoods.h"
#include "../TQualityMap.h"

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
	double estimate(TPopulationLikehoodLocus & storage, double epsilonF, TGlfConverter & glfConverter);
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
	double LL_atMAP;
	double f_CI_lower, f_CI_upper;
	double* f_initialGrid;
	double* f_grid;
	double* LL_initialGrid;
	double* Likelihood_grid;


	double guessInitialAlleleFrequency(TPopulationLikehoodLocus & storage, TGlfConverter & glfConverter);
	double calcLL(TPopulationLikehoodLocus & storage, THardyWeinbergGenotypeProbabilities & pGenotype, TGlfConverter & glfConverter);
	void fillInitialGrid(TPopulationLikehoodLocus & storage, TGlfConverter & glfConverter);
	void estimateMAP(TPopulationLikehoodLocus & storage, TGlfConverter & glfConverter);
	void estimateCredibleIntervals(TPopulationLikehoodLocus & storage, TGlfConverter & glfConverter);

public:
	TAlleleFreqEstimatorBayes(TParameters & Parameters, TLog* logfile, TRandomGenerator* RandomGenerator);
	~TAlleleFreqEstimatorBayes();
	double estimate(TPopulationLikehoodLocus & storage, TGlfConverter & glfConverter);

	double credibleIntervalUsed(){ return credibleInterval; };
	double MAP(){ return f_MAP; };
	double lowerCredibleInterval(){ return f_CI_lower; };
	double upperCredibleInterval(){ return f_CI_upper; };
};

//------------------------------------------------
//TAlleleFreqEstimator
//------------------------------------------------
class TAlleleFreqEstimator{
	// about vcf-file
	std::string vcfFilename;
	bool vcfRead;
	TLog* logfile;

	// data on individuals
	TPopulationSamples samples;
	TGlfConverter glfConverter;

	//data on loci
	long _numLoci;

public:
	TAlleleFreqEstimator(TParameters & Parameters, TLog* logfile);
	void estimateAlleleFreq(TParameters & Parameters, TRandomGenerator* randomGenerator);
};



#endif /* POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_ */
