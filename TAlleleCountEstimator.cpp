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
	alleleFrequencyLikelihoods_h = new double[numAlleleCounts];
};

TSiteAlleleFrequencyLikelihoods::~TSiteAlleleFrequencyLikelihoods(){
	delete[] alleleFrequencyLikelihoods_h;
};

/*
void TSiteAlleleFrequencyLikelihoods::fill(unsigned short* phred){
	//Calculating allele frequency likelihoods according to Nielsen et al. (2012) PLoS One, page 3
	//initialize
	alleleFrequencyLikelihoods_h[0] = phred[0];
	alleleFrequencyLikelihoods_h[1] = phred[1];
	alleleFrequencyLikelihoods_h[2] = phred[2];

	for(int j=3; j<numAlleleCounts; j++)
		alleleFrequencyLikelihoods_h[j] = 0.0;

	//Recursion
	for(int d=2; d<numInd_k; ++d){
		for(int j=2*(d+1); j>1; j--){

			alleleFrequencyLikelihoods_h[j]
		}
	}

	//Termination

};
*/

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
