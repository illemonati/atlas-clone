/*
 * TInbreedingEstimator.cpp
 *
 *  Created on: Dec 6, 2018
 *      Author: linkv
 */

#include "TInbreedingEstimator.h"

//---------------------------
// allele frequencies p
//---------------------------
TAlleleFreq::TAlleleFreq(){
	numLoci = -1;
	proposalWidths = NULL;
	sumIterations = NULL;
	sumOfSquaresIterations = NULL;
}


TAlleleFreq::TAlleleFreq(std::vector<double> & P, float initialProposalWidth){
	alleleFreq = P;
	numLoci = alleleFreq.size();
	proposalWidths = new float[numLoci];
	sumIterations = new double[numLoci];
	sumOfSquaresIterations = new double[numLoci];
	//init proposal widths and check range of p
	for(int l=0; l<numLoci; ++l){
		proposalWidths[l] = initialProposalWidth;
		if(alleleFreq[l] == 0)
			//TODO: is it ok to approximate like this?
			//small number
			alleleFreq[l] = 0.000000000001;
		else if(alleleFreq[l] == 1)
			alleleFreq[l] = 0.99999999999;
	}
}

void TAlleleFreq::setToValue(double fixedValue){
	for(int l=0; l<numLoci; ++l){
		alleleFreq[l] = fixedValue;
	}
}

void TAlleleFreq::adjustProposalWidthAfterBurnin(int* numAcceptedP, int numUpdates){
	//adjust proposal with for p
	for(long l=0; l<numLoci; ++l){
		double newProposalWidth = proposalWidths[l];
		newProposalWidth *=  (double) numAcceptedP[l] / ((double) numUpdates * 3.0);

		if(newProposalWidth / proposalWidths[l] < 0.1)
			newProposalWidth = 0.1 * proposalWidths[l];
		else if(proposalWidths[l] / newProposalWidth < 0.1)
			newProposalWidth = 10 * proposalWidths[l];

		proposalWidths[l] = newProposalWidth;

		if(proposalWidths[l] <= 0)
			throw "Proposal width for allele frequency at locus " + toString(l) + " is not larger than 0!";

	}
//	numUpdates = 0;
}

void TAlleleFreq::update(long & index, double & value){
	if(value < 0)
		throw "updating allele freq at locus " + toString(index) + " to negative value!";
	if(value == 0)
		//TODO: is it ok to approximate like this?
		//small number
		value = 0.000000000001;
	else if(value == 1)
		value = 0.99999999999;
	alleleFreq[index] = value;
}

double TAlleleFreq::proposeNew(long & locusNum, TRandomGenerator & randomGenerator){
	double newP = alleleFreq[locusNum] + randomGenerator.getRand() * proposalWidths[locusNum] - proposalWidths[locusNum] / 2.0;
	//is proposed p outside range [0,1]
	while(newP > 1 || newP < 0){
		//mirror
		if(newP < 0)
			newP = -newP;
		if(newP > 1)
			newP = 2 - newP;
	}
	return newP;
}

//TAlleleFreq::updateAll(TPopulationLikelihoods & likelihoods, TAlphaOrBeta & alpha, TAlphaOrBeta & beta, TRandomGenerator & randomGenerator){
//	for(int l=0; l<numLoci; ++l){
//		//propose new p
//		double newP = p[l] + randomGenerator.getRand() * proposalWidths[l] - proposalWidths[l] / 2.0;
//		while(newP > 1 || newP < 0){
//			if(newP < 0)
//				newP = -newP;
//			if(newP > 1)
//				newP = 2 - newP;
//		}
//
//		double sumOverInds = 0.0;
//
//		//calculate hastings ratio
//		for(likelihoods.begin(); !likelihoods.end(); likelihoods.next()){
//			uint8_t* data = likelihoods.curData();
//			for(int s=0; s<likelihoods.curSampleSize(); ++s){
//				int index = 3*s;
//				//calculate and add ratio for each genotype
//				for(int g=0; g<3; ++g){
//					sumOverInds += log((data[index + g] * PGenoGivenFAndP(g, F, newP)) /  (data[index + g] * PGenoGivenFAndP(g, F, p[l])));
//				}
//			}
//		}
//
//		double logH = (alpha.getNaturalScaleValue() - 1) * (log(newP) - log(p[l]))
//					+ (beta.getNaturalScaleValue() - 1) * (log(1 - newP) - log(1 - p[l]))
//					+ logLikelihoodAllInds(newP, F, alpha, beta)
//					- logLikelihoodAllInds(p[l], F, alpha, beta);
//
//		//accept?
//		if(log(randomGenerator.getRand()) < logH){
//			p[l] = newP;
//			return true;
//		} else
//			return false;
//	}
//}

//---------------------------
// alphaOrBeta
//---------------------------

TAlphaOrBeta::TAlphaOrBeta(){
	variableName = "";
	proposalWidth = -1.0;
	_alphaOrBeta = -1.0;
	_logAlphaOrBeta = -1.0;
}

TAlphaOrBeta::TAlphaOrBeta(std::string VariableName, double & ProposalWidth){
	variableName = VariableName;
	proposalWidth = ProposalWidth;
	_alphaOrBeta = -1.0;
	_logAlphaOrBeta = -1.0;
}

void TAlphaOrBeta::update(double & newLogValue, double & newNaturalScaleValue){
	if(newNaturalScaleValue < 0)
		throw "updating " + variableName + " to negative value!";
	_logAlphaOrBeta = newLogValue;
	_alphaOrBeta = newNaturalScaleValue;
}

void TAlphaOrBeta::adjustProposalWidthAfterBurnin(int numAccepted, int numUpdates){
	//adjust proposal for alpha / beta
	if(numAccepted == 0)
		throw "numAccepted = 0";
	double newProposalWidth = proposalWidth;
	newProposalWidth *= (double) numAccepted / 2.0 / (double) numUpdates * 3.0;

	if(newProposalWidth / proposalWidth < 0.1)
		newProposalWidth = 0.1 * proposalWidth;
	else if(proposalWidth / newProposalWidth < 0.1)
		newProposalWidth = 10 * proposalWidth;

	proposalWidth = newProposalWidth;

	if(proposalWidth <= 0)
		throw "Proposal width for parameter " + variableName + " is not larger than 0!";
}

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

	widthProposalKernelP = Parameters.getParameterDoubleWithDefault("widthProposalKernelP", 0.05);
	logfile->list("Will use a proposal kernel of width " + toString(widthProposalKernelP) + " for updates of log(alpha) and log(beta)");

	numBurnins = Parameters.getParameterIntWithDefault("numBurnins", 1);
	burninLength = Parameters.getParameterIntWithDefault("burninLength", 1000);
	logfile->list("Will run " + toString(numBurnins) + " burnin(s) of length " + toString(burninLength) + " and adjust the proposal kernel widths after each before starting the full MCMC");

	thinning = Parameters.getParameterIntWithDefault("thinning", 1);
	if(thinning < 1 || thinning > numIterations)
		throw "Thinning must be > 1 and < number iterations!";
	if(thinning > 1){
		if(thinning == 2)
			logfile->list("Will print every second iterations to the output file (thinning = 2)");
		else if(thinning == 3)
			logfile->list("Will print every third iterations to the output file (thinning = 3)");
		else
			logfile->list("Will print every " + toString(thinning) + "th iterations to the output file (thinning = " + toString(thinning) + ")");
	}

	//read data
	likelihoods.doSaveAlleleFrequencies();
	if(Parameters.parameterExists("trueAlleleFreq")){
		likelihoods.doSaveTrueAlleleFrequencies();
		logfile->list("Will save true allele frequencies in populationLikelihoods object");
	}
	likelihoods.readData(Parameters, logfile);
	numLoci = likelihoods.getNumLoci();

	//initialize parameters
	initParams(randomGenerator, Parameters);
	numAcceptedF = 0;
	numAcceptedP = new int[numLoci];
	for(unsigned int l=0; l<numLoci; ++l)
		numAcceptedP[l] = 0;
	numAcceptedAlpha = 0;
	numAcceptedBeta = 0;

	//TODO: get defaultOutName from vcf file in TPopulationLikelihoods
	std::string defaultOutName = "default";
	outname = Parameters.getParameterStringWithDefault("outname", defaultOutName);
}

void TInbreedingEstimator::initializeAlphaAndBeta(){
	// estimate alpha_f and beta_f by method of moments
	double mean = 0.0;
	double sumXSquare = 0.0;
	double alphaF;
	double betaF;

	for(int l = 0; l < numLoci; l++){
		mean += p[l];
		sumXSquare += p[l] * p[l];
	}
	mean /= (double) numLoci;
	sumXSquare /= (double) numLoci;
	double var = sumXSquare - mean  * mean;

	//now estimate alpha and beta
	double tmp = ((mean * (1.0 - mean)) / var ) - 1.0;
	alphaF = mean * tmp;
	if(alphaF < 0.0)
		alphaF = 0.01;
	double logAlpha = log(alphaF);
	alpha.update(logAlpha, alphaF);
	betaF = (1.0 - mean) * tmp;
	if(betaF < 0.0)
		betaF = 0.01;
	double logBeta = log(betaF);
	beta.update(logBeta, betaF);
	logfile->list("Initialized alpha to " + toString(alphaF) + " and beta to " + toString(betaF));
}

void TInbreedingEstimator::initParams(TRandomGenerator & randomGenerator, TParameters & parameters){
	//F
//	F = 0.0;
	double tmp = randomGenerator.getRand();
	if(tmp <= pi)
		F = 0;
	else
		F = randomGenerator.getRand();
	F = exp(randomGenerator.getRand(0, 1));
	logfile->list("initialized F to " + toString(F));

	//p
	std::vector<double> tmp2;
	if(parameters.parameterExists("trueAlleleFreq")){
		tmp2 = likelihoods.donateTrueAlleleFrequencies();
		logfile->list("initializing allele frequencies to true values read from trueAlleleFreq file");
	} else {
		tmp2 = likelihoods.donateAlleleFrequencies();
		logfile->list("initializing allele frequencies to values estimated from sample genotype likelihoods");
	}
	p = TAlleleFreq(tmp2, widthProposalKernelP);
	p.setToValue(0.2);

	if(p.numLoci != numLoci)
		throw "Did not receive one allele frequency per locus! Number of loci=" + toString(numLoci) + " and number of allele frequencies " + toString(p.numLoci);
	logfile->list("initialized allele frequencies for " + toString(p.numLoci) + " loci");

	//alpha
	double widthProposalKernelLogAlpha = parameters.getParameterDoubleWithDefault("widthProposalKernelLogAlpha", 0.35);
	alpha = TAlphaOrBeta("alpha", widthProposalKernelLogAlpha);

	//beta
	double widthProposalKernelLogBeta = parameters.getParameterDoubleWithDefault("widthProposalKernelLogBeta", 0.35);
	beta = TAlphaOrBeta("beta", widthProposalKernelLogBeta);

	initializeAlphaAndBeta();

	logfile->list("Will use a proposal kernel of width " + toString(widthProposalKernelLogAlpha) + " for updates of log(alpha)");
	logfile->list("Will use a proposal kernel of width " + toString(widthProposalKernelLogBeta) + " for updates of log(beta)");
}

void TInbreedingEstimator::printTrajectory(gz::ogzstream & tracefile){
	tracefile << F << "\t" << alpha.getNaturalScaleValue() << "\t" << beta.getNaturalScaleValue() << "\n";
}

bool TInbreedingEstimator::updateF(){
	std::cout << "updating F" << std::endl;
	return true;
}

bool TInbreedingEstimator::updateP(uint8_t* data, long & locusNum, int curSampleSize, TAlphaOrBeta & alpha, TAlphaOrBeta & beta){
	//propose new p
	double newP = p.proposeNew(locusNum, randomGenerator);

	//calculate hastings ratio
	double sumOverInds = 0.0;
	for(int s=0; s<likelihoods.curSampleSize(); ++s){
		int index = 3*s;
		//calculate and add ratio for each genotype
		for(int g=0; g<3; ++g){
			std::cout << "data for g:" << g << " is " << (unsigned) data[index + g] << " exp: " << exp(data[index + g]) << std::endl;
			sumOverInds += log((qualMap.phredToError(data[index + g]) * probGenoGivenFAndP(g, F, newP)) /  (qualMap.phredToError(data[index + g]) * probGenoGivenFAndP(g, F, p[locusNum])));
		}
	}

	double logH = (alpha.getNaturalScaleValue() - 1) * (log(newP) - log(p[locusNum]))
				+ (beta.getNaturalScaleValue() - 1) * (log(1 - newP) - log(1 - p[locusNum]))
				+ logLikelihoodAllInds(data, likelihoods.curSampleSize(), newP, F, alpha, beta)
				- logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[locusNum], F, alpha, beta);

	if(logLikelihoodAllInds(data, likelihoods.curSampleSize(), newP, F, alpha, beta) - logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[locusNum], F, alpha, beta) != sumOverInds){
		std::cout << logLikelihoodAllInds(data, likelihoods.curSampleSize(), newP, F, alpha, beta) - logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[locusNum], F, alpha, beta) << std::endl;
		std::cout << sumOverInds << std::endl;
		std::cout << logH << std::endl;
		std::cout << "alpha.getNaturalScaleValue() " << alpha.getNaturalScaleValue() << std::endl;
		std::cout << "beta.getNaturalScaleValue() " << beta.getNaturalScaleValue() << std::endl;
		std::cout << "newP " << newP << std::endl;
		std::cout << "log(1 - p[locusNum]) " << log(1 - p[locusNum]) << std::endl;
		int g = 2;
		std::cout << "probGenoGivenFAndP(g, F, newP) " << probGenoGivenFAndP(g, F, newP) << std::endl;

		throw "two values are not equal!";
	}

	//accept?
	if(log(randomGenerator.getRand()) < logH){
		p.update(locusNum, newP);
		return true;
	} else
		return false;
}

bool TInbreedingEstimator::updateAlphaOrBeta(TAlphaOrBeta & alphaOrBetaToUpdate, TAlphaOrBeta & alphaOrBetaOther){
	// compute sum
	double sumP = 0;
	for (unsigned int l = 0; l < numLoci; l++){
		sumP += log(p[l]);
	}

	// propose new log(value) (uniform proposal kernel)
	double newLogAlphaOrBeta = alphaOrBetaToUpdate.getLogValue() + randomGenerator.getRand() * alphaOrBetaToUpdate.proposalWidth - alphaOrBetaToUpdate.proposalWidth / 2.0;
	double newAlphaOrBeta = exp(newLogAlphaOrBeta);
	if(newAlphaOrBeta < 0)
		throw "proposed new log(alpha) that on natural scale is negative!";

	// compute log hastings ratio
	double logH = numLoci * (randomGenerator.gammaln(newAlphaOrBeta + alphaOrBetaOther.getNaturalScaleValue())
			+ randomGenerator.gammaln(alphaOrBetaToUpdate.getNaturalScaleValue())
			- randomGenerator.gammaln(alphaOrBetaToUpdate.getNaturalScaleValue() + alphaOrBetaOther.getNaturalScaleValue())
			- randomGenerator.gammaln(newAlphaOrBeta))
			+ (newAlphaOrBeta - alphaOrBetaToUpdate.getNaturalScaleValue()) * sumP;

//	if(alphaOrBetaToUpdate.variableName == "alpha"){
//		std::cout << std::endl;
//		std::cout << "logAlpha " << alphaOrBetaToUpdate.getLogValue() << std::endl;
//		std::cout << "alphaOrBetaToUpdate.proposalWidth " << alphaOrBetaToUpdate.proposalWidth << std::endl;
//		std::cout << "newLogAlpha: " << newLogAlphaOrBeta << std::endl;
//		std::cout << "logH " << logH << std::endl;
//	}
	// accept or reject
	if(log(randomGenerator.getRand()) < logH){
		alphaOrBetaToUpdate.update(newLogAlphaOrBeta, newAlphaOrBeta);
		return true;
	}
	else
		return false;
}

double TInbreedingEstimator::probGenoGivenFAndP(int & genotype, double & F, double & p){
	double probGenoGivenFAndP;
	if(genotype == 0)
		probGenoGivenFAndP = (1.0 - F)*p*p + F*p;
	else if(genotype == 1)
		probGenoGivenFAndP = 2.0*p*(1.0 - p)*(1.0 - F);
	else if(genotype == 2)
		probGenoGivenFAndP = (1.0 - F)*(1.0 - p)*(1.0 - p) + F*(1.0 - p);
	else
		throw "unknown genotype '" + toString(genotype) +"'!";
	if(probGenoGivenFAndP > 1.0 || probGenoGivenFAndP < 0.0)
		throw "probGenoGivenFAndP is not between 0 and 1! It is " + toString(probGenoGivenFAndP) + ". genotype is " + toString(genotype) + " p is " + toString(p) + " and F is " + toString(F);
	return probGenoGivenFAndP;
}

double TInbreedingEstimator::logLikelihoodAllInds(uint8_t* data, int curSampleSize, double & thisP, double & thisF, TAlphaOrBeta & alpha, TAlphaOrBeta & beta){
	//sum over all individuals of log sum_g P(d|g)P(g|p,F)
	double sumOverInds = 0.0;

	double PGeno[3];
	PGeno[0] = (1.0 - thisF)*(1.0 - thisP)*(1.0 - thisP) + thisF*(1.0 - thisP);
	PGeno[1] = 2.0*thisP*(1.0 - thisP)*(1.0 - thisF);
	PGeno[2] = (1.0 - thisF)*thisP*thisP + thisF*thisP;

	for(int s=0; s<curSampleSize; ++s){
		int index = 3*s;
		//calculate and add ratio for each genotype


		//TODO: find way to exclude samples where all g.t. likelihoods = 0
		//check if we should include sample
//		double sumOfPhredToError = 0;
//		for(int g=0; g<3; ++g)
//			sumOfPhredToError += data[index + g];
//		if(sumOfPhredToError == 0.0)
//			continue;

		//sample has data -> add it to likelihood

		double integrationOverGeno = qualMap.phredToError(data[index]) * PGeno[0];
		integrationOverGeno += qualMap.phredToError(data[index + 1]) * PGeno[1];
		integrationOverGeno += qualMap.phredToError(data[index + 2]) * PGeno[2];

		//check if likelihood of sample is a probability
		if(integrationOverGeno <= 0)
			throw "likelihood is not larger than 0!";
//		if(integrationOverGeno > 1){
//			for(int g=0; g<3; ++g)
//				std::cout << "g " << g << " qualMap.phredToError(data[index + g] " << qualMap.phredToError(data[index + g]) << std::endl;
//			throw "integration over genotypes is larger than 1!";
//			-> likelihoods don't have to sum to 1!
//		}

		sumOverInds += log(integrationOverGeno);
	}

	if(sumOverInds != sumOverInds )
		throw "sumOverInds is not a number!";
//	if(sumOverInds > 0)
//		throw "logLikelihoodAllInds is larger than zero! logLikelihoodAllInds=" + toString(sumOverInds);
	return sumOverInds;
}

double TInbreedingEstimator::logProbPGivenAlphaBeta(){
	//P(p|alpha, beta)
	double posteriorProbability = 0.0;
	//calc sum of log(f) and log(1-f)
	double sumLogFreq = 0.0;
	double sumLogOneMinusFreq = 0.0;

	for(unsigned long l=0; l<numLoci; ++l){
		sumLogFreq += log(p[l]);
		sumLogOneMinusFreq += log(1.0 - p[l]);
	}

	//add beta density
	posteriorProbability += (alpha.getNaturalScaleValue()-1.0) * sumLogFreq + (beta.getNaturalScaleValue() - 1.0) * sumLogOneMinusFreq;
	posteriorProbability += numLoci * (randomGenerator.gammaln(alpha.getNaturalScaleValue() + beta.getNaturalScaleValue()) - randomGenerator.gammaln(alpha.getNaturalScaleValue()) - randomGenerator.gammaln(beta.getNaturalScaleValue()));

//	if(posteriorProbability > 0)
//		throw "logProbPGivenAlphaBeta is larger than 1!";
	return posteriorProbability;
}
void TInbreedingEstimator::writeLikelihoodForDebuggingF(TParameters & params){
	//open output file
	std::string tracefile = outname + "_logLikelihoodF.txt.gz";
	logfile->list("Will write likelihood chain for grid of F to file '" + tracefile + "'.");
	gz::ogzstream outP(tracefile.c_str());
	if(!outP)
		throw "Failed to open file '" + tracefile + "' for writing!";
	std::cout << "first three true alleleFreq: " << p[0] << " " <<  p[1] << " " << p[2] << std::endl;

	//initialize
	double trueAlpha = 0.5;
	double trueLogAlpha = log(trueAlpha);
	alpha.update(trueLogAlpha, trueAlpha);
	double trueBeta = 0.5;
	double trueLogBeta = log(trueBeta);
	beta.update(trueLogBeta, trueBeta);

	//make grid
	int numProbs = 100;
	float step = 1.0 / (float) numProbs;
	std::cout << "step size "  << step << std::endl;

	//calculate ll
	for(int i=0; i<numProbs; ++i){
		double logLikelihood = 0.0;
		double thisF = 0.0;
		thisF = thisF + (double) i*step;
		if(thisF > 1)
			break;

		std::cout << "thisF " << thisF << std::endl;
		long l = 0;
		for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
			uint8_t* data = likelihoods.curData();
			logLikelihood += logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], thisF, alpha, beta);
		}

		std::cout << "logLikelihood before adding posterior " << logLikelihood << std::endl;

		//add P(p|alpha,beta)
//		double tmp =logProbPGivenAlphaBeta();
//		std::cout << "likelihood individuals " << logLikelihood << "posterior prob " << tmp << std::endl;
		logLikelihood += logProbPGivenAlphaBeta();
		outP << thisF << "\t" << logLikelihood << "\n";
		std::cout << logLikelihood << std::endl;
	}
	logfile->done();
}

void TInbreedingEstimator::writeLikelihoodForDebuggingAlleleFreq(TParameters & params){
	//which locus to update
	long index = 0;

	//open output file
	std::string tracefile = outname + "_logLikelihoodP[" + toString(index) + "].txt.gz";
	logfile->list("Will write likelihood chain for grid of p[" + toString(index) + "] to file '" + tracefile + "'.");
	logfile->list("The true value for p[" + toString(index) + "] is: " + toString(p[index]));
	gz::ogzstream outP(tracefile.c_str());
	if(!outP)
		throw "Failed to open file '" + tracefile + "' for writing!";
	std::cout << "first three true alleleFreq: " << p[0] << " " <<  p[1] << " " << p[2] << std::endl;

	//initialize
	double thisF = 0.0;
	double trueAlpha = 0.5;
	double trueLogAlpha = log(trueAlpha);
	alpha.update(trueLogAlpha, trueAlpha);
	double trueBeta = 0.5;
	double trueLogBeta = log(trueBeta);
	beta.update(trueLogBeta, trueBeta);

	//make grid
	int numProbs = 100;
	float step = 1.0 / (float) numProbs;
	double pValue = 0.0;

	//calculate ll
	for(int i=0; i<numProbs; ++i){
		double newPValue = pValue + (double) i*step;
		std::cout << "calculating likelihood for p[0]=" << newPValue << std::endl;
		p.update(index, newPValue);
		double logLikelihood = 0.0;
		long l = 0;
		for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
			uint8_t* data = likelihoods.curData();
			logLikelihood += logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], thisF, alpha, beta);
			if(logLikelihood > 0)
				throw "likelihood is larger than 1!";
		}

		std::cout << "logLikelihood before adding posterior " << logLikelihood << std::endl;

		//add P(p|alpha,beta)
//		double tmp =logProbPGivenAlphaBeta();
//		std::cout << "likelihood individuals " << logLikelihood << "posterior prob " << tmp << std::endl;
		logLikelihood += logProbPGivenAlphaBeta();
		outP << newPValue << "\t" << logLikelihood << "\n";
		std::cout << logLikelihood << std::endl;
	}
	logfile->done();
}

void TInbreedingEstimator::writeLikelihoodForDebuggingAlpha(TParameters & params){
	//open output file
	std::string tracefile = outname + "_logLikelihoodAlpha.txt.gz";
	logfile->list("Will write likelihoods for grid of alpha to '" + tracefile + "'.");
	gz::ogzstream outP(tracefile.c_str());
	if(!outP)
		throw "Failed to open file '" + tracefile + "' for writing!";

	//initialize
	double thisF = 0.0;
	double alphaValue = 0.0;
	int numProbs = 100;
	float step = 2.0 / (float) numProbs;
	double trueBeta = 0.5;
	double trueLogBeta = log(trueBeta);
	beta.update(trueLogBeta, trueBeta);
	std::cout << "first three true alleleFreq: " << p[0] << " " <<  p[1] << " " << p[2] << std::endl;

	//calculate ll
	for(int i=0; i<numProbs; ++i){
		double newAlphaValue = alphaValue + (double) i*step;
		std::cout << "calculating likelihood for alpha=" << newAlphaValue << std::endl;
		if(newAlphaValue == 0.0)
			//a small number
			newAlphaValue = 0.000000001;
		double newAlphaValueLog = log(newAlphaValue);
		alpha.update(newAlphaValueLog, newAlphaValue);
		double logLikelihood = 0.0;
		long l = 0;
		for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
			uint8_t* data = likelihoods.curData();
			logLikelihood += logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], thisF, alpha, beta);
			if(logLikelihood > 0)
				throw "likelihood is larger than 1!";
		}

		//add P(p|alpha,beta)
//		double tmp =logProbPGivenAlphaBeta();
//		std::cout << "likelihood individuals " << logLikelihood << "posterior prob " << tmp << std::endl;
		logLikelihood += logProbPGivenAlphaBeta();
		outP << alpha.getNaturalScaleValue() << "\t" << logLikelihood << "\n";
//		std::cout << logLikelihood << std::endl;
	}
	logfile->done();
}

void TInbreedingEstimator::oneMCMCIteration(int iterationNum){
	//update params
//	numAcceptedF += updateF();
	//locus index

	long l = 0;
	for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
		uint8_t* data = likelihoods.curData();
		numAcceptedP[l] += updateP(data, l, likelihoods.curSampleSize(), alpha, beta);
	}
	//alpha
	numAcceptedAlpha += updateAlphaOrBeta(alpha, beta);
	//beta
	numAcceptedBeta += updateAlphaOrBeta(beta, alpha);
}

void TInbreedingEstimator::printAcceptanceRates(int numIterations){
	logfile->conclude("acceptance rate for F is " + toString((double) numAcceptedF / numIterations));
	logfile->conclude("acceptance rate for alpha is " + toString((double) numAcceptedAlpha / numIterations));
	logfile->conclude("acceptance rate for beta is " + toString((double) numAcceptedBeta / numIterations));
	logfile->conclude("acceptance rate for first locus is " + toString((double) numAcceptedP[0] / numIterations));
}

void TInbreedingEstimator::resetAcceptanceRates(){
	numAcceptedF = 0;
	numAcceptedAlpha = 0;
	numAcceptedBeta = 0;
	for(unsigned int l=0; l<numLoci; ++l)
		numAcceptedP[l] = 0;
}

void TInbreedingEstimator::adjustProposalWidths(){
//	F.adjustProposalWidth(numAcceptedF, numIterations);
	p.adjustProposalWidthAfterBurnin(numAcceptedP, burninLength);
	alpha.adjustProposalWidthAfterBurnin(numAcceptedAlpha, burninLength);
	beta.adjustProposalWidthAfterBurnin(numAcceptedBeta, burninLength);
}

void TInbreedingEstimator::writeParameterEstimatesOfIteration(gz::ogzstream & out){
	out << F << "\t" << alpha.getNaturalScaleValue() << "\t" << alpha.getLogValue() << "\t" << beta.getNaturalScaleValue() << beta.getLogValue() << "\n";
}

void TInbreedingEstimator::runEstimation(){
	logfile->startIndent("Running MCMC to estimate inbreeding coefficient F and the allele frequency distribution:");

	//open MCMC output file
	std::string tracefile = outname + "_inbreedingMCMC.txt.gz";
	logfile->list("Will write MCMC chain to file '" + tracefile + "'.");
	gz::ogzstream out(tracefile.c_str());
	if(!out)
		throw "Failed to open file '" + tracefile + "' for writing!";

	tracefile = outname + "_inbreedingMCMC_pSelectedLoci.txt.gz";
	logfile->list("Will write MCMC chain for selected loci to file '" + tracefile + "'.");
	gz::ogzstream outP(tracefile.c_str());
	if(!outP)
		throw "Failed to open file '" + tracefile + "' for writing!";


	//---------------------------
	//run Burnin(s)
	//---------------------------
	for(int b=0; b<numBurnins; ++b){
		std::string reportString = "Running a burnin of " + toString(burninLength) + " iterations (";
		logfile->listFlush(reportString + "0%) ...");
		int oldProg=0;
		for(long i=0; i<burninLength; ++i){

			oneMCMCIteration(i);

			//report
			int prog = floor((float) i / (float) burninLength * 100);
			if(prog > oldProg){
				oldProg = prog;
				logfile->listOverFlush(reportString + toString(prog) + "%) ...");
			}
		}
		logfile->overList("Running a burnin of " + toString(burninLength) + " iterations ... done!  ");

		//print acceptance rates
		printAcceptanceRates(burninLength);

		//adjust things after burnin
		logfile->listFlush("Adjusting proposal widths ...");
		adjustProposalWidths();
		logfile->done();

		resetAcceptanceRates();
	}

	//---------------------------
	//run MCMC
	//---------------------------

	//write header
	out << "F\talpha\talphaLog\tbeta\tbetaLog\n";

	//write initial parameter estimates
	writeParameterEstimatesOfIteration(out);

	//run MCMC
	int oldProg = 0;
	std::string progressString = "Running MCMC chain of length " + toString(numIterations) + " ...";
	logfile->listFlush(progressString + "(0%)");
	for(int i=1; i<numIterations; ++i){

		oneMCMCIteration(i);

		//print to file
		if(i % thinning == 0){
			writeParameterEstimatesOfIteration(out);
			outP << p[10] << "\t" << p[15] << "\n";
		}

		// print progress
		int prog = (double) i / (double) numIterations * 100.0;
		if(prog > oldProg){
			oldProg = prog;
			logfile->listOverFlush(progressString + " (" + toString(oldProg) + "%)");
		}
	}
	logfile->overList(progressString + " done!   ");
	printAcceptanceRates(numIterations);

	//clean up
	logfile->endIndent();
	out.close();
}
