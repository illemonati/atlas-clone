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
	TQualityMap phredToGTLMap;

public:
	TAlleleFreqEstimatorHardyWeinberg();
	double estimate(TPopulationLikehoodStorage & storage, double epsilonF);
};

//------------------------------------------------
//TAlleleFreqEstimatorBayes
//------------------------------------------------
class TAlleleFreqEstimatorBayes{
private:
	TRandomGenerator* randomGenerator;
	int burninLength;
	int numBurnins;
	std::vector<double> mcmc;

	double alpha, beta;
	double oneMinusAlpha, oneMinusBeta;

	double f, f_lower, f_upper;
	int lowerIndex, upperIndex;

	TQualityMap phredToGTLMap;

	double guessInitialAlleleFrequency(TPopulationLikehoodStorage & storage);
	double calcLL(TPopulationLikehoodStorage & storage, THardyWeinbergGenotypeProbabilities & pGenotype);
	int makeMCMCUpdate(TPopulationLikehoodStorage & storage, double & oldLL, const double & prop, THardyWeinbergGenotypeProbabilities* pGenotype, int & old);

public:
	TAlleleFreqEstimatorBayes(TParameters & Parameters, TLog* logfile, TRandomGenerator* RandomGenerator);
	double estimate(TPopulationLikehoodStorage & storage);

	double posteriorMean(){ return f; };
	double lowerCredibleInterval(){ return f_lower; };
	double upperCredibleInterval(){ return f_upper; };
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

	//data on loci
	long _numLoci;

public:
	TAlleleFreqEstimator(TParameters & Parameters, TLog* logfile);
	void estimateAlleleFreq(TParameters & Parameters, TRandomGenerator* randomGenerator);
};



#endif /* POPULATIONTOOLS_TALLELEFREQUENCYESTIMATOR_H_ */
