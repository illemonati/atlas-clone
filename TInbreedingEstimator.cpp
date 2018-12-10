/*
 * TInbreedingEstimator.cpp
 *
 *  Created on: Dec 6, 2018
 *      Author: linkv
 */

#include "TInbreedingEstimator.h"

TInbreedingEstimator::TInbreedingEstimator(TParameters & Parameters, TLog* Logfile){
	logfile = Logfile;

	//algorithm params
	numIterations = Parameters.getParameterIntWithDefault("numIter", 1000);
	logfile->list("Stopping MCMC after " + toString(numIterations) + " interations");

	pi = Parameters.getParameterDoubleWithDefault("pi", 0.1);
	logfile->list("Mixture parameter pi set to " + toString(pi));

	//read data
	likelihoods.doSaveAlleleFrequencies();
	likelihoods.readData(Parameters, logfile);

	//initialize
	initParams(randomGenerator);
	//TODO: get defaultOutName from vcf file in TPopulationLikelihoods
	std::string defaultOutName = "default";
	outname = Parameters.getParameterStringWithDefault("outname", defaultOutName);

}

void TInbreedingEstimator::initializeAlphaBeta(){
	// estimate alpha_f and beta_f by method of moments
	double mean = 0.0;
	double sumXSquare = 0.0;

	for(unsigned int i=0; i<p.size(); ++i){
		mean += p[i];
		sumXSquare += p[i] * p[i];
	}

	mean /= (double) p.size();
	sumXSquare /= (double) p.size();
	double var = sumXSquare - mean  * mean;

	//now estimate alpha and beta
	double tmp = ((mean * (1.0 - mean)) / var ) - 1.0;
	alpha = mean * tmp;
	if(alpha < 0.0)
		alpha = 0.01;
	logAlpha = log(alpha);
	beta = (1.0 - mean) * tmp;
	if(beta < 0.0)
		beta = 0.01;
	logBeta= log(beta);
	logfile->list("Initialized alpha to " + toString(alpha) + " and beta to " + toString(beta));
}

void TInbreedingEstimator::initParams(TRandomGenerator & randomGenerator){
	double tmp = randomGenerator.getRand();
	if(tmp <= pi)
		F = 0;
	else
		F = randomGenerator.getRand();
	F = exp(randomGenerator.getRand(0, 1));
	p = likelihoods.donateAlleleFrequencies();  //new double[likelihoods.getNumLoci()];
	initializeAlphaBeta();
}

void TInbreedingEstimator::printTrajectory(gz::ogzstream & tracefile){
	tracefile << F << "\t" << alpha << "\t" << beta << "\n";
}

void TInbreedingEstimator::updateF(){
	std::cout << "updating F" << std::endl;
}

void TInbreedingEstimator::updateP(long l){
}

void TInbreedingEstimator::updateAlphaOrBeta(){
//
//	// compute sum
//	double sumF = 0;
//	for (int l = 0; l < numLoci; l++){
//		sumF += log(logisticLookup->approxLogistic(loci->mus[l]));
//	}
//
//	// propose new log(alphaF) (uniform proposal kernel)
//	double newLogAlphaF = loci->logAlphaF + randomGenerator->getRand() * widthProposalKernelLogAlphaF
//			- widthProposalKernelLogAlphaF / 2.0;
//	double newAlphaF = exp(newLogAlphaF);
//
//	// compute log hastings ratio
//	double logH = numLoci * (randomGenerator->gammaln(newAlphaF + loci->betaF)
//			- randomGenerator->gammaln(loci->alphaF + loci->betaF)
//			+ randomGenerator->gammaln(loci->alphaF)
//			- randomGenerator->gammaln(newAlphaF))
//			+ (newAlphaF - loci->alphaF) * sumF;
//
//	// accept or reject
//	if (log(randomGenerator->getRand()) < logH){
//		loci->setAlphaF(newAlphaF); // set new alphaF and logAlphafF
//		loci->setLogAlphaF(newLogAlphaF);
//		return true;
//	}
//	else
//		return false;
}


void TInbreedingEstimator::runEstimation(){
	logfile->startIndent("Running MCMC to estimate inbreeding coefficient F and the allele frequency distribution:");

	//open MCMC output file
	std::string tracefile = outname + "_inbreedingMCMC.txt.gz";
	logfile->list("Will write MCMC chain to file '" + tracefile + "'.");
	gz::ogzstream out(tracefile.c_str());
	if(!out)
		throw "Failed to open file '" + tracefile + "' for writing!";

	//write header
	out << "F\talpha\tbeta\n";

	//run MCMC
	int oldProg = 0;
	for(int i=0; i<numIterations; ++i){
		std::string progressString = "Running MCMC chain of length " + toString(numIterations) + " ...";
		logfile->listFlush(progressString + "(0%)");

		//update params
		updateF();
		//locus index
		long l = 0;
		for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
			updateP(l);
		}
		//alpha
		updateAlphaOrBeta();
		//beta
		updateAlphaOrBeta();

		// print progress
		int prog = (double) i / (double) numIterations * 100.0;
		if(prog > oldProg){
			oldProg = prog;
			logfile->listOverFlush(progressString + " (" + toString(oldProg) + "%)");
		}
	}

	//clean up
	logfile->endIndent();
	out.close();
}
