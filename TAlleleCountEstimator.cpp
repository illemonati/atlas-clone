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

	//fill inverse of choose
	log_choose_2k_j = new double[numAlleleCounts];
	for(int j=0; j<numAlleleCounts; ++j)
		log_choose_2k_j[j] = math::chooseLog(2*numInd_k, j);
};

TSiteAlleleFrequencyLikelihoods::~TSiteAlleleFrequencyLikelihoods(){
	delete[] alleleFrequencyLikelihoods_h;
};


/*
double angsd::addProtect2(double a,double b){
  //function does: log(exp(a)+exp(b)) while protecting for underflow
  double maxVal;// = std::max(a,b));
  if(a>b)
    maxVal=a;

  else
    maxVal=b;
  double sumVal = exp(a-maxVal)+exp(b-maxVal);
  return log(sumVal) + maxVal;
}


double angsd::addProtect3(double a,double b, double c){
  //function does: log(exp(a)+exp(b)+exp(c)) while protecting for underflow
  double maxVal;// = std::max(a,std::max(b,c));
  if(a>b&&a>c)
    maxVal=a;
  else if(b>c)
    maxVal=b;
  else
    maxVal=c;
  double sumVal = exp(a-maxVal)+exp(b-maxVal)+exp(c-maxVal);
  return log(sumVal) + maxVal;
}

void TSiteAlleleFrequencyLikelihoods::fill(unsigned short* phred){
	//Calculating allele frequency likelihoods according to Nielsen et al. (2012) PLoS One, page 3
	//initialize
	alleleFrequencyLikelihoods_h[0] = phredToGTLMap[ phred[0] ];
	alleleFrequencyLikelihoods_h[1] = phredToGTLMap[ phred[0] ];
	alleleFrequencyLikelihoods_h[2] = phredToGTLMap[ phred[0] ];

	for(int j=3; j<numAlleleCounts; j++)
		alleleFrequencyLikelihoods_h[j] = 0.0;

	//Recursion
	for(int d=2; d<numInd_k; ++d){
		int s = 3*d;
		for(int j=2*(d+1); j>1; j--){
			alleleFrequencyLikelihoods_h[j] = phredToGTLMap[ phred[s + 2] ] * alleleFrequencyLikelihoods_h[j-2]
										    + phredToGTLMap[ phred[s + 1] ] * alleleFrequencyLikelihoods_h[j-1]
											+ phredToGTLMap[ phred[s] ] * alleleFrequencyLikelihoods_h[j];
		}

		//special case for j=1,0
		alleleFrequencyLikelihoods_h[1] = phredToGTLMap[ phred[s + 1] ] * alleleFrequencyLikelihoods_h[0]
									    + phredToGTLMap[ phred[s] ] * alleleFrequencyLikelihoods_h[1];
		alleleFrequencyLikelihoods_h[0] = phredToGTLMap[ phred[s] ] * alleleFrequencyLikelihoods_h[0];
	}

	//Termination
	for(int j=0; j<numAlleleCounts; j++)
		alleleFrequencyLikelihoods_h[j] = log(alleleFrequencyLikelihoods_h[j]) - log_choose_2k_j[j];


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
