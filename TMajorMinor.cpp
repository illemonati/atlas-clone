/*
 * TMajorMinor.cpp
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#include "TMajorMinor.h"

TMajorMinor::TMajorMinor(TLog* Logfile){
	logfile = Logfile;
};



void TMajorMinor::fillInitialEstimateOfGenotypeFrequencies(double* genoFreq, const int & numSamples, double* genotypeLikelihoods){
	//calculate by using MLE genotype for each individual
	genoFreq[0] = 0.0; genoFreq[1] = 0.0; genoFreq[2] = 0.0;
	for(int i = 0 ; i < 3 * numSamples; i += 3){
		if(genotypeLikelihoods[i + 1] > genotypeLikelihoods[i]){
			if(genotypeLikelihoods[i + 2] > genotypeLikelihoods[i + 1]) genoFreq[2] += 1.0;
			else genoFreq[1] += 1.0;
		} else {
			if(genotypeLikelihoods[i + 2] > genotypeLikelihoods[i]) genoFreq[2] += 1.0;
			else genoFreq[0] += 1.0;
		}
		/*
		for(int g = 0; g<3; ++g){
			genoFreq[g] += genotypeLikelihoods[i][g];
		}
		*/
	}

	double sum = 0.0;
	for(int g = 0; g<3; ++g){
		genoFreq[g] /= (double) numSamples;
		if(genoFreq[g] <= 0.0) genoFreq[g] = 0.01;
		if(genoFreq[g] >= 1.0) genoFreq[g] = 0.99;
		sum += genoFreq[g];
	}
	for(int g = 0; g<3; ++g){
		genoFreq[g] /= sum;
	}
};

void TMajorMinor::estimateGenotypeFrequencies(double* genotypeFrequencies, const int & numSamples, double* genotypeLikelihoods, const double & epsilonF){
	//prepare variables
	double sum;
	int i, g;
	double weightsNull[3];
	double genoFreq_old[3];

	//estimate initial frequencies from MLEs
	fillInitialEstimateOfGenotypeFrequencies(genotypeFrequencies, numSamples, genotypeLikelihoods);

	//run EM for max 1000 steps
	for(int s=0; s<1000; ++s){
		//set genofreq and calc P(g|f)
		for(g=0; g<3; ++g){
			genoFreq_old[g] = genotypeFrequencies[g];
			genotypeFrequencies[g] = 0.0;
		}

		//estimate new genotype frequencies
		for(i = 0; i < 3 * numSamples; i += 3){
			sum = 0.0;
			for(g = 0; g < 3; ++g){
				weightsNull[g] = genotypeLikelihoods[i + g] * genoFreq_old[g];
				sum += weightsNull[g];
			}
			genotypeFrequencies[0] += weightsNull[0] / sum;
			genotypeFrequencies[2] += weightsNull[2] / sum;
		}

		genotypeFrequencies[0] /= (double) numSamples;
		genotypeFrequencies[2] /= (double) numSamples;
		genotypeFrequencies[1] = 1.0 - genotypeFrequencies[0] - genotypeFrequencies[2];

		//check if we stop
		if(fabs(genotypeFrequencies[0] - genoFreq_old[0]) < epsilonF && fabs(genotypeFrequencies[2] - genoFreq_old[2]) < epsilonF) break;
	}
};

void TMajorMinor::calculateLogLikelihood(const int & numSamples, double* genotypeLikelihoods, double* genotypeFrequencies){

}

void TMajorMinor::fillLikelihoods(const int & numSamples, uint8_t** phred, double* genotypeLikelihoods, int aa, int ab, int bb){
	for(int i=0; i<numSamples; ++i){
		int ind = 3*i;
		genotypeLikelihoods[ind] = qualiMap[ phred[i][aa] ];
		genotypeLikelihoods[ind + 1] = qualiMap[ phred[i][ab] ];
		genotypeLikelihoods[ind + 2] = qualiMap[ phred[i][bb] ];
	}
}

void TMajorMinor::getMLEOfMajorMinor(const int & numSamples, uint8_t** phred, double* genotypeLikelihoods, const double & epsilonF){
	//prepare variables
	double genotypeFrequencies[3];
	double LL[6];

	//run across all six allele combinations: A/C, A/G, A/T, C/G, C/T and G/T
	// A/C
	fillLikelihoods(numSamples, phred, genotypeLikelihoods, 0, 1, 4);
	estimateGenotypeFrequencies(genotypeFrequencies, numSamples, genotypeLikelihoods, epsilonF);





};


void TMajorMinor::estimateMajorMinor(TParameters & params){
	//open GLF files
	TGlfMultiReader glfReader(params, logfile);
	glfReader.setAllActive();

	//TODO: think about filters

	//variables
	double* genotypeLikelihoods = new double[3*glfReader.numSamples()];

	while(glfReader.readNext()){
		//DEBUG: print
		glfReader.print();



		//do major minor




	}
};
