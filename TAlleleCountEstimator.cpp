/*
 * TAlleleCountEstimator.cpp
 *
 *  Created on: Dec 9, 2018
 *      Author: phaentu
 */


#include "TAlleleCountEstimator.h"



//-------------------------------------------------
// TSiteAlleleFrequencyLikelihoods
//-------------------------------------------------
TSiteAlleleFrequencyLikelihoods::TSiteAlleleFrequencyLikelihoods(int numIndividuals){
	numInd_k = numIndividuals;
	numAlleleCounts = 2*numInd_k + 1;
	log_alleleFrequencyLikelihoods_h = new double[numAlleleCounts];

	//fill inverse of choose
	log_choose_2k_j = new double[numAlleleCounts];
	for(int j=0; j<numAlleleCounts; ++j)
		log_choose_2k_j[j] = chooseLog(2*numInd_k, j);
};

TSiteAlleleFrequencyLikelihoods::~TSiteAlleleFrequencyLikelihoods(){
	delete[] log_alleleFrequencyLikelihoods_h;
};


double TSiteAlleleFrequencyLikelihoods::protectedSumInLog(double a,double b){
  //returns log(exp(a)+exp(b)) while protecting for underflow, inspired by ANGSD
  double maxVal;
  if(a>b) maxVal = a;
  else maxVal = b;
  double sumVal = exp(a-maxVal)+exp(b-maxVal);
  return log(sumVal) + maxVal;
};

double TSiteAlleleFrequencyLikelihoods::protectedSumInLog(double a, double b, double c){
  //returns log(exp(a)+exp(b)+exp(c)) while protecting for underflow, inspired by ANGSD
  double maxVal;
  if(a > b && a > c) maxVal = a;
  else if(b > c) maxVal = b;
  else maxVal = c;
  double sumVal = exp(a-maxVal)+exp(b-maxVal)+exp(c-maxVal);
  return log(sumVal) + maxVal;
};

void TSiteAlleleFrequencyLikelihoods::fill(unsigned short* phred){
	//Calculating allele frequency likelihoods according to Nielsen et al. (2012) PLoS One, page 3
	//initialize
	log_alleleFrequencyLikelihoods_h[0] = phred[0];
	log_alleleFrequencyLikelihoods_h[1] = phred[0];
	log_alleleFrequencyLikelihoods_h[2] = phred[0];

	for(int j=3; j<numAlleleCounts; j++)
		log_alleleFrequencyLikelihoods_h[j] = 0.0;

	//Recursion
	for(int d=2; d<numInd_k; ++d){
		int s = 3*d;
		for(int j=2*(d+1); j>1; j--){
			log_alleleFrequencyLikelihoods_h[j] = protectedSumInLog(
											phred[s + 2] + log_alleleFrequencyLikelihoods_h[j-2],
										    phred[s + 1] + log_alleleFrequencyLikelihoods_h[j-1],
											phred[s] + log_alleleFrequencyLikelihoods_h[j] );
		}

		//special case for j=1,0
		log_alleleFrequencyLikelihoods_h[1] = protectedSumInLog(
											phred[s + 1] + log_alleleFrequencyLikelihoods_h[0],
											phred[s] + log_alleleFrequencyLikelihoods_h[1] );
		log_alleleFrequencyLikelihoods_h[0] = phred[s] + log_alleleFrequencyLikelihoods_h[0];
	}

	//Termination
	for(int j=0; j<numAlleleCounts; j++)
		log_alleleFrequencyLikelihoods_h[j] = log(log_alleleFrequencyLikelihoods_h[j]) - log_choose_2k_j[j];


};

//-------------------------------------------------
// TAlleleCountEstimator
//-------------------------------------------------
TAlleleCountEstimator::TAlleleCountEstimator(TParameters & params, TLog* Logfile){
	logfile = Logfile;

	//initialize random generator
	//TODO: do the random generator initialization in the task switcher?
	logfile->listFlush("Initializing random generator ...");
	if(params.parameterExists("fixedSeed")){
		randomGenerator = new TRandomGenerator(params.getParameterLong("fixedSeed"), true);
	} else if(params.parameterExists("addToSeed")){
		randomGenerator = new TRandomGenerator(params.getParameterLong("addToSeed"), false);
	} else randomGenerator = new TRandomGenerator();
	logfile->write(" done with seed " + toString(randomGenerator->usedSeed) + "!");
};

TAlleleCountEstimator::~TAlleleCountEstimator(){
	delete randomGenerator;
};

void TAlleleCountEstimator::estimateAlleleCounts(TParameters & params){

	//first do a test: read Pop likelihoods and print them again.
	TPopulationLikelihoods likelihoods(params, logfile);

	likelihoods.print();



};
