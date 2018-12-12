/*
 * TInbreedingEstimator.cpp
 *
 *  Created on: Dec 6, 2018
 *      Author: linkv
 */

#include "TInbreedingEstimator.h"

//---------------------------
// alphaOrBeta
//---------------------------

TAlphaOrBeta::TAlphaOrBeta(std::string VariableName){
	variableName = VariableName;
	_alphaOrBeta = -1.0;
	_logAlphaOrBeta = -1.0;
}

bool TAlphaOrBeta::initialize(std::vector<double> & p, TLog* logfile){
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
	_alphaOrBeta = mean * tmp;
	if(_alphaOrBeta < 0.0)
		_alphaOrBeta = 0.01;
	_logAlphaOrBeta = log(_alphaOrBeta);

	logfile->list("Initialized " + variableName + " to " + toString(_alphaOrBeta));
	return true;
}

void TAlphaOrBeta::update(double & newLogValue, double & newNaturalScaleValue){
	_logAlphaOrBeta = newLogValue;
	_alphaOrBeta = newNaturalScaleValue;
}

//void alphaOrBeta::updateWithNaturalScaleValue(double & newValue){
//	_alphaOrBeta = newValue;
//	_logAlphaOrBeta = log(newValue);
//}

double TAlphaOrBeta::getLogValue(){
	return _logAlphaOrBeta;
}

double TAlphaOrBeta::getNaturalScaleValue(){
	return _alphaOrBeta;
}
//---------------------------
// TInbreedingEstimator
//---------------------------

TInbreedingEstimator::TInbreedingEstimator(TParameters & Parameters, TLog* Logfile){
	logfile = Logfile;

	//algorithm params
	numIterations = Parameters.getParameterIntWithDefault("numIter", 1000);
	logfile->list("Stopping MCMC after " + toString(numIterations) + " interations");

	pi = Parameters.getParameterDoubleWithDefault("pi", 0.1);
	logfile->list("Mixture parameter pi set to " + toString(pi));

	widthProposalKernelLogAlphaOrBeta = Parameters.getParameterDoubleWithDefault("widthProposalKernelLogAlphaAndBeta", 0.35);
	logfile->list("Will use a proposal kernel of width " + toString(widthProposalKernelLogAlphaOrBeta) + " for updates of log(alpha) and log(beta)");


	//read data
	likelihoods.doSaveAlleleFrequencies();
	likelihoods.readData(Parameters, logfile);
	numLoci = likelihoods.getNumLoci();

	//initialize
	initParams(randomGenerator);
	//TODO: get defaultOutName from vcf file in TPopulationLikelihoods
	std::string defaultOutName = "default";
	outname = Parameters.getParameterStringWithDefault("outname", defaultOutName);

}

void TInbreedingEstimator::initParams(TRandomGenerator & randomGenerator){
	//F
	double tmp = randomGenerator.getRand();
	if(tmp <= pi)
		F = 0;
	else
		F = randomGenerator.getRand();
	F = exp(randomGenerator.getRand(0, 1));

	//p
	p = likelihoods.donateAlleleFrequencies();  //new double[likelihoods.getNumLoci()];
	if(p.size() != numLoci)
		throw "Did not receive one allele frequency per locus! Number of loci=" + toString(numLoci) + " and number of allele frequencies " + toString(p.size());

	//alpha
	alpha = TAlphaOrBeta("alpha");
	if(!alpha.initialize(p, logfile))
		throw "failed to initialize alpha!";

	//beta
	beta = TAlphaOrBeta("beta");
	if(!beta.initialize(p, logfile))
		throw "failed to initialize beta!";
}

void TInbreedingEstimator::printTrajectory(gz::ogzstream & tracefile){
	tracefile << F << "\t" << alpha.getNaturalScaleValue() << "\t" << beta.getNaturalScaleValue() << "\n";
}

void TInbreedingEstimator::updateF(){
	std::cout << "updating F" << std::endl;
}

void TInbreedingEstimator::updateP(long l, TAlphaOrBeta & alpha, TAlphaOrBeta & beta){
	//propose new p
	double newP;
	double sumOverInds = 0.0;

	//calculate hastings ratio
	for(likelihoods.begin(); !likelihoods.end(); likelihoods.next()){
		uint8_t* data = likelihoods.curData();
		for(int s=0; s<likelihoods.curSampleSize(); ++s){
			int index = 3*s;
			//calculate and add ratio for each genotype
			for(int g=0; g<3; ++g){
				sumOverInds += log((data[index + g] * PGenoGivenFAndP(g, F, newP)) /  (data[index + g] * PGenoGivenFAndP(g, F, p[l])));
			}
		}
	}

	double h = sumOverInds + (alpha.getNaturalScaleValue() - 1) * (log(newP) - log(p[l]))
			+ (beta.getNaturalScaleValue() - 1) * (log(1 - newP) - log(1 - p[l]));
	//accept?
}

bool TInbreedingEstimator::updateAlphaOrBeta(TAlphaOrBeta & alphaOrBetaToUpdate, TAlphaOrBeta & alphaOrBetaOther){
	// compute sum
	double sumP = 0;
	for (unsigned int l = 0; l < numLoci; l++){
		sumP += log(p[l]);
	}

	// propose new log(value) (uniform proposal kernel)
	double newLogAlphaOrBeta = alphaOrBetaToUpdate.getLogValue() + randomGenerator.getRand() * widthProposalKernelLogAlphaOrBeta - widthProposalKernelLogAlphaOrBeta / 2.0;
	double newAlphaOrBeta = exp(newLogAlphaOrBeta);

	// compute log hastings ratio
	double logH = numLoci * (randomGenerator.gammaln(newAlphaOrBeta + alphaOrBetaOther.getNaturalScaleValue())
			+ randomGenerator.gammaln(alphaOrBetaToUpdate.getNaturalScaleValue())
			- randomGenerator.gammaln(alphaOrBetaToUpdate.getNaturalScaleValue() + alphaOrBetaOther.getNaturalScaleValue())
			- randomGenerator.gammaln(newAlphaOrBeta))
			+ (newAlphaOrBeta - alphaOrBetaToUpdate.getNaturalScaleValue()) * sumP;

	// accept or reject
	if (log(randomGenerator.getRand()) < logH){
		alphaOrBetaToUpdate.update(newLogAlphaOrBeta, newAlphaOrBeta);
		return true;
	}
	else
		return false;
}

double TInbreedingEstimator::PGenoGivenFAndP(int & genotype, double & F, double & p){
	if(genotype == 0)
		return (1 - F)*p*p + F;
	else if(genotype == 1)
		return 2*p*(1 - p)*(1 - F);
	else if(genotype == 2)
		return (1 - F)*(1 - p)*(1 - p) + F*(1 - p);
	else
		throw "unknown genotype '" + toString(genotype) +"'!";
}

void TInbreedingEstimator::oneMCMCIteration(){
	//update params
	updateF();
	//locus index
	long l = 0;
	for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
		updateP(l, alpha, beta);
	}
	//alpha
	updateAlphaOrBeta(alpha, beta);
	//beta
	updateAlphaOrBeta(beta, alpha);

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

		oneMCMCIteration();

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
