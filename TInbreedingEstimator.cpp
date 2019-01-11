/*
 * TInbreedingEstimator.cpp
 *
 *  Created on: Dec 6, 2018
 *      Author: linkv
 */

#include "TInbreedingEstimator.h"

//---------------------------
// F
//---------------------------
TInbreedingF::TInbreedingF(){
	_F = -1.0;
	_probMovingToModelNoF = -1.0;
	_sdProposal = -1.0;
	_inModelWithF = false;
	_lambda = -1.0;
	_posteriorProbModelWithF = -1;
}

TInbreedingF::TInbreedingF(double F, float & ProbMovingToModelNoF, double & SdProposal, bool InModelWithF, double Lambda){
	_F = F;
	_probMovingToModelNoF = ProbMovingToModelNoF;
	_sdProposal = SdProposal;
	_inModelWithF = InModelWithF;
	_lambda = Lambda;
	_posteriorProbModelWithF = 0;
}

void TInbreedingF::adjustProposalWidthAfterBurnin(int numAcceptedF, int numUpdates, TLog* logfile){
//	if(numAcceptedF == 0)
//		throw "numAccepted = 0";
	double newProposalWidth = _sdProposal;
	newProposalWidth *=  ((double) numAcceptedF / (double) numUpdates) * 3.0;

	if(newProposalWidth / _sdProposal < 0.1)
		newProposalWidth = 0.1 * _sdProposal;
	else if(_sdProposal / newProposalWidth < 0.1)
		newProposalWidth = 10 * _sdProposal;

	if(newProposalWidth > 1.0)
		newProposalWidth = 0.1;

	if(newProposalWidth < 0.00001)
		newProposalWidth = 0.1;

	_sdProposal = newProposalWidth;

	if(!(_sdProposal > 0))
		throw "standard deviation of F proposal kernel is not larger than 0!";

}

double TInbreedingF::proposeNew(TRandomGenerator & randomGenerator){
    double newF = _F + randomGenerator.getNormalRandom(0, _sdProposal);
    while(newF > 1 || newF < 0){
        //mirror
        if(newF < 0)
            newF = -newF;
        if(newF > 1)
            newF = 2 - newF;
    }
    if(newF < 0)
    	throw "proposing F smaller than zero!";
    if(newF > 1)
    	throw "proposing F larger than one!";
    return newF;
}

void TInbreedingF::update(double value, bool inModelWithF){
	_F = value;
	_inModelWithF = inModelWithF;
	_posteriorProbModelWithF += inModelWithF;
}

double TInbreedingF::logPDFExp(){
	return log(_lambda) - _lambda * _F;
}

double TInbreedingF::PDFExp(){
	return _lambda * exp(-_lambda* _F);
}

float TInbreedingF::probMovingToModelNoF(){
	return _probMovingToModelNoF;
}

double TInbreedingF::F(){
	return _F;
}

bool TInbreedingF::inModelWithF(){
	return _inModelWithF;
}

double TInbreedingF::lambda(){
	return _lambda;
}

int TInbreedingF::posteriorProbModelWithF(){
	return _posteriorProbModelWithF;
}

//---------------------------
// allele frequencies p
//---------------------------
TAlleleFreq::TAlleleFreq(){
	numLoci = -1;
	proposalWidths = NULL;
	sumIterations = NULL;
	sumOfSquaresIterations = NULL;
	minAlleleFreq = -1.0;
}


TAlleleFreq::TAlleleFreq(std::vector<double> & P, float initialProposalWidthFactor, int numSamples){
	alleleFreq = P;
	numLoci = alleleFreq.size();
	proposalWidths = new float[numLoci];
	sumIterations = new double[numLoci];
	sumOfSquaresIterations = new double[numLoci];
	minAlleleFreq = 1.0 / (double) numSamples / 10.0;
	//init proposal widths and check range of p
	for(int l=0; l<numLoci; ++l){
		sumIterations[l] = 0.0;
		sumOfSquaresIterations[l] = 0.0;
		if(alleleFreq[l] < minAlleleFreq)
			//TODO: is it ok to approximate like this?
			//small number
			alleleFreq[l] = minAlleleFreq;
		else if(alleleFreq[l] > 1.0 - minAlleleFreq)
			alleleFreq[l] = 1.0 - minAlleleFreq;
		proposalWidths[l] = initialProposalWidthFactor * alleleFreq[l];

	}
}

void TAlleleFreq::setSumsToZero(){
	for(int l=0; l<numLoci; ++l){
		sumIterations[l] = 0.0;
		sumOfSquaresIterations[l] = 0.0;
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
		if(numAcceptedP[l] == 0){
			std::cout << "!!!!!!!!! acceptance rate is zero for locus " << l << std::endl;
//			double minFreq = 1.0 / (double) numSamples;
//			alleleFreq[l] = minFreq;
			proposalWidths[l] *= 0.1;
			continue;
		}

		double newProposalWidth = proposalWidths[l];
		newProposalWidth *=  (double) numAcceptedP[l] / (double) numUpdates * 3.0;

		if(newProposalWidth / proposalWidths[l] < 0.1)
			newProposalWidth = 0.1 * proposalWidths[l];
		else if(proposalWidths[l] / newProposalWidth < 0.1)
			newProposalWidth = 10.0 * proposalWidths[l];

		if(newProposalWidth > 0.5)
			newProposalWidth = 0.5;

		if(l == 0)
			std::cout << "old proposalWidths[0] " << proposalWidths[l] << " newProposalWidth " << newProposalWidth << std::endl;


		proposalWidths[l] = newProposalWidth;


		if(proposalWidths[l] <= 0)
			throw "Proposal width for allele frequency at locus " + toString(l) + " is not larger than 0!";

	}
//	numUpdates = 0;
}

void TAlleleFreq::update(long & index, double value){
	if(value < 0)
		throw "updating allele freq at locus " + toString(index) + " to negative value!";
	if(alleleFreq[index] < minAlleleFreq)
		//TODO: is it ok to approximate like this?
		//small number
		alleleFreq[index] = minAlleleFreq;
	else if(alleleFreq[index] > 1.0 - minAlleleFreq)
		alleleFreq[index] = 1.0 - minAlleleFreq;
	alleleFreq[index] = value;
	sumIterations[index] += value;
	sumOfSquaresIterations[index] += value*value;
}

double TAlleleFreq::proposeNew(long & locusNum, TRandomGenerator & randomGenerator){
//	newAlleleFreq[G][l] = fabs(alleleFreq[G][l] + (randomGenerator.getRand()-0.5) * proposalWidthFrequencies[G][l]);

	double newP = alleleFreq[locusNum] + randomGenerator.getRand() * proposalWidths[locusNum] - proposalWidths[locusNum] / 2.0;

	//is proposed p outside range [0,1]
	if(newP < minAlleleFreq){
//			std::cout << "mirroring p at 0: " << newP << " to " << -newP << std::endl;
		newP = 2.0*minAlleleFreq - newP;
	} else if(newP > 1.0 - minAlleleFreq){
//			std::cout << "mirroring p at 1: " << newP << " to " << 2 - newP << std::endl;
			newP = 2.0 * (1.0 - minAlleleFreq) - newP;
	}
	return newP;
}

double TAlleleFreq::getPosteriorMean(unsigned long & index, int numUpdates){
	return sumIterations[index] / (double) numUpdates;
}

double TAlleleFreq::getPosteriorVariance(unsigned long & index, int numUpdates){
	return sumOfSquaresIterations[index] / (double) numUpdates - getPosteriorMean(index, numUpdates);
}

//---------------------------
// gamma
//---------------------------

TGamma::TGamma(){
	proposalWidth = -1.0;
	_gamma = -1.0;
	_logGamma = -1.0;
}

TGamma::TGamma(double & ProposalWidth){
	proposalWidth = ProposalWidth;
	_gamma = -1.0;
	_logGamma = -1.0;
}

void TGamma::update(const double & newLogValue, const double & newNaturalScaleValue){
	if(newNaturalScaleValue < 0)
		throw "updating gamma to negative value!";
	_logGamma = newLogValue;
	_gamma = newNaturalScaleValue;
}

void TGamma::adjustProposalWidthAfterBurnin(int numAccepted, int numUpdates){
	//adjust proposal for alpha / beta
//	if(numAccepted == 0)
//		throw "numAccepted = 0";
	double newProposalWidth = proposalWidth;
	newProposalWidth *= (double) numAccepted / (double) numUpdates * 3.0;

	if(newProposalWidth / proposalWidth < 0.1)
		newProposalWidth = 0.1 * proposalWidth;
	else if(proposalWidth / newProposalWidth < 0.1)
		newProposalWidth = 10 * proposalWidth;
	if(newProposalWidth > 1.0)
		newProposalWidth = 1.0;

	proposalWidth = newProposalWidth;

	if(proposalWidth <= 0)
		throw "Proposal width for gamma is not larger than 0!";
}

double TGamma::getLogValue(){
	return _logGamma;
}

double TGamma::getNaturalScaleValue(){
	return _gamma;
}
//---------------------------
// TInbreedingEstimator
//---------------------------

TInbreedingEstimator::TInbreedingEstimator(TParameters & Parameters, TLog* Logfile){
	logfile = Logfile;

	//algorithm params
	numIterations = Parameters.getParameterIntWithDefault("numIter", 1000);
	logfile->list("Stopping MCMC after " + toString(numIterations) + " interations");

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
	numAcceptedP = new int[numLoci];
	resetAcceptanceRates();

	//output name
	std::string vcfFileName = likelihoods.getVCFName();
	vcfFileName = extractBeforeLast(vcfFileName, ".vcf");
	outname = Parameters.getParameterStringWithDefault("out", vcfFileName);
}

void TInbreedingEstimator::initializeGamma(){
	// estimate alpha_f and beta_f by method of moments
	double mean = 0.0;
	double sumXSquare = 0.0;
	double alphaF;
//	double betaF;

	for(unsigned int l = 0; l < numLoci; l++){
		mean += p[l];
		sumXSquare += p[l] * p[l];
	}
	mean /= (double) numLoci;
	sumXSquare /= (double) numLoci;
	double var = sumXSquare - mean  * mean;

	//now estimate gamma -> just assume it's same as alpha
	double tmp = ((mean * (1.0 - mean)) / var ) - 1.0;
	alphaF = mean * tmp;
	if(alphaF < 0.0)
		alphaF = 0.01;
	double logAlpha = log(alphaF);
	Gamma.update(logAlpha, alphaF);
//	betaF = (1.0 - mean) * tmp;
//	if(betaF < 0.0)
//		betaF = 0.01;
//	double logBeta = log(betaF);
//	beta.update(logBeta, betaF);



	double trueVal = 0.2;
	double logTrueVal = log(0.2);
	Gamma.update(logTrueVal, trueVal);
//	beta.update(logTrueVal, trueVal);

	logfile->list("Initialized gamma to " + toString(Gamma.getNaturalScaleValue()));
}

void TInbreedingEstimator::initParams(TRandomGenerator & randomGenerator, TParameters & parameters){
	//F
	sdF = parameters.getParameterDoubleWithDefault("sdF", 0.2);
	logfile->list("Standard deviation of proposal kernel for F is set to " + toString(sdF));

	float probMovingToModelNoF = parameters.getParameterDoubleWithDefault("probMovingToModelNoF", 0.1);
	logfile->list("Will propose move to model without F with probability " + toString(probMovingToModelNoF));

	double lambda = parameters.getParameterDoubleWithDefault("lambda", 1.0);
	logfile->list("Lambda of exponential distribution used for the proposal of new F after move to model with F is set to " + toString(lambda));

	bool startInModelWithF = parameters.parameterExists("startInModelWithF");
	if(startInModelWithF){
		F = TInbreedingF(randomGenerator.getExponentialRandom(lambda), probMovingToModelNoF, sdF, startInModelWithF, lambda);
		F.update(randomGenerator.getExponentialRandomTruncated(lambda, 0, 1), true);
		F.update(0.2, true);
		logfile->list("initialized F to " + toString(F.F()) + " in model " + toString(F.inModelWithF()));
	} else {
		F = TInbreedingF(0, probMovingToModelNoF, sdF, startInModelWithF, lambda);
		logfile->list("initialized F to " + toString(F.F()) + " in model " + toString(F.inModelWithF()));
	}

	double widthProposalKernelP = parameters.getParameterDoubleWithDefault("widthProposalKernelPFactor", 1);
	logfile->list("Will use a proposal kernel of width " + toString(widthProposalKernelP) + " for updates of log(alpha) and log(beta)");

	//p
	std::vector<double> tmp2;
	if(parameters.parameterExists("trueAlleleFreq")){
		tmp2 = likelihoods.donateTrueAlleleFrequencies();
		trueAlleleFreq = tmp2;
		trueAlleleFreqProvided = true;
		logfile->list("initializing allele frequencies to true values read from trueAlleleFreq file");
	} else {
		tmp2 = likelihoods.donateAlleleFrequencies();
		logfile->list("initializing allele frequencies to values estimated from sample genotype likelihoods");
		trueAlleleFreqProvided = false;
	}
	p = TAlleleFreq(tmp2, widthProposalKernelP, likelihoods.getNumIndividuals());
	if(parameters.parameterExists("trueAlleleFreq"))
			std::cout << "first 5 true p: " << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << "," << p[4] << std::endl;

//	p.setToValue(0.2);

	if(p.numLoci != numLoci)
		throw "Did not receive one allele frequency per locus! Number of loci=" + toString(numLoci) + " and number of allele frequencies " + toString(p.numLoci);
	logfile->list("initialized allele frequencies for " + toString(p.numLoci) + " loci");

	//gamma
	double widthProposalKernelGamma = parameters.getParameterDoubleWithDefault("widthProposalKernelGamma", 0.35);
	Gamma = TGamma(widthProposalKernelGamma);

	initializeGamma();

	logfile->list("Will use a proposal kernel of width " + toString(widthProposalKernelGamma) + " for updates of log(alpha)");
	logfile->list("Will use a proposal kernel of width " + toString(widthProposalKernelGamma) + " for updates of log(beta)");
}

bool TInbreedingEstimator::updateF(){
	if(F.inModelWithF()){
		//try move to model without F
		if(randomGenerator.getRand() < F.probMovingToModelNoF()){
			double logH = F.logPDFExp() - log(F.probMovingToModelNoF());

			long l = 0;
			for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
				uint8_t* data = likelihoods.curData();
				logH += logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], 0)
				- logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], F.PDFExp());
			}

			//accept?
			if(log(randomGenerator.getRand()) < logH){
				F.update(0, false);
				return true;
			} else
				return false;

		} else {
			//propose new F
			double newF = F.proposeNew(randomGenerator);

			long l = 0;
			double logH = 0.0;
			for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
				uint8_t* data = likelihoods.curData();
				logH += logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], newF)
				- logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], F.F());
			}

			//accept?
			if(log(randomGenerator.getRand()) < logH){
				F.update(newF, true);
				return true;
			} else
				return false;
		}
	}

	//moving to model with F
	else {
		double logH = log(F.probMovingToModelNoF()) - F.logPDFExp();
		double newF = randomGenerator.getExponentialRandomTruncated(F.lambda(), 0.0, 1.0);

		long l = 0;
		for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
			uint8_t* data = likelihoods.curData();
			logH += logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], newF)
			- logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], 0);
		}

		//accept?
		if(log(randomGenerator.getRand()) < logH){
			F.update(newF, true);
			return true;
		} else
			return false;
	}
}

bool TInbreedingEstimator::updateP(uint8_t* data, long & locusNum, int curSampleSize, TGamma & Gamma){
	//propose new p
	double newP = p.proposeNew(locusNum, randomGenerator);

//	Gamma.update(log(0.2), 0.2);

	double logH = (Gamma.getNaturalScaleValue() - 1.0) * (log(newP) + log(1.0 - newP) - log(p[locusNum]) - log(1.0 - p[locusNum])) // - log(1.0 - p[locusNum]));													)
				+ logLikelihoodAllInds(data, likelihoods.curSampleSize(), newP, F.F())
				- logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[locusNum], F.F());


//	double logHold = (Gamma.getNaturalScaleValue() - 1.0) * (log(newP) - log(p[locusNum]))
//				+ (Gamma.getNaturalScaleValue() - 1.0) * (log(1 - newP) - log(1 - p[locusNum]))
//				+ logLikelihoodAllInds(data, likelihoods.curSampleSize(), newP, F.F())
//				- logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[locusNum], F.F());
//
//	std::cout << "logH " << logH << " logHold " << logHold << std::endl;

	//accept?
	if(log(randomGenerator.getRand()) < logH){
		p.update(locusNum, newP);
		return true;
	} else {
		p.update(locusNum, p[locusNum]);
		return false;
	}
}

bool TInbreedingEstimator::updateGamma(){
	// propose new log(value) (uniform proposal kernel)
	double newLogGamma = Gamma.getLogValue() + randomGenerator.getRand() * Gamma.proposalWidth - Gamma.proposalWidth / 2.0;

//	Gamma.update(log(0.5),0.5);
//	double newLogGamma = 0.5;

	double newGamma = exp(newLogGamma);
	if(newGamma < 0)
		throw "proposed new log(alpha) that on natural scale is negative!";

	//calc sum of log(f) and log(1-f)
	double sumLogFreq = 0.0;
	double sumLogOneMinusFreq = 0.0;

	for(unsigned long l=0; l<numLoci; ++l){
		sumLogFreq += log(p[l]);
		sumLogOneMinusFreq += log(1.0 - p[l]);
	}

	double tmpH = numLoci * (randomGenerator.gammaln(newGamma + Gamma.getNaturalScaleValue())
			+ randomGenerator.gammaln(Gamma.getNaturalScaleValue())
			- randomGenerator.gammaln(Gamma.getNaturalScaleValue() + Gamma.getNaturalScaleValue())
			- randomGenerator.gammaln(newGamma))
			+ (newGamma - Gamma.getNaturalScaleValue()) * sumLogFreq;

//	std::cout << "tmpH " << tmpH << std::endl;


	// compute log hastings ratio
	double diffGammas = newGamma - Gamma.getNaturalScaleValue();
	double logH = (double) numLoci *
			(randomGenerator.gammaln(2.0*newGamma)
			+ 2.0*randomGenerator.gammaln(Gamma.getNaturalScaleValue())
			- randomGenerator.gammaln(2.0*Gamma.getNaturalScaleValue())
			- 2.0*randomGenerator.gammaln(newGamma)
			)
			+ diffGammas * sumLogFreq
			+ diffGammas * sumLogOneMinusFreq;

//	std::cout << "logH " << logH << std::endl;

	// accept or reject
	double tmp = log(randomGenerator.getRand());
	if(tmp < logH){
		Gamma.update(newLogGamma, newGamma);
		return true;
	}
	else
		return false;
}

double TInbreedingEstimator::logLikelihoodAllInds(uint8_t* data, int curSampleSize, double thisP, double thisF){
	//sum over all individuals of log sum_g P(d|g)P(g|p,F)
	double sumOverInds = 0.0;

	double PGeno[3];
	PGeno[0] = (1.0 - thisF)*(1.0 - thisP)*(1.0 - thisP) + thisF*(1.0 - thisP);
	PGeno[1] = 2.0*thisP*(1.0 - thisP)*(1.0 - thisF);
	PGeno[2] = (1.0 - thisF)*thisP*thisP + thisF*thisP;

	for(int s=0; s<curSampleSize; ++s){
		int index = 3*s;
		//calculate and add ratio for each genotype
		double integrationOverGeno = qualMap.phredToError(data[index]) * PGeno[0];
		integrationOverGeno += qualMap.phredToError(data[index + 1]) * PGeno[1];
		integrationOverGeno += qualMap.phredToError(data[index + 2]) * PGeno[2];

		//check if likelihood of sample is a probability
		if(integrationOverGeno <= 0)
			throw "likelihood is not larger than 0!";

		sumOverInds += log(integrationOverGeno);
	}

	if(sumOverInds != sumOverInds )
		throw "sumOverInds is not a number!";
	return sumOverInds;
}

double TInbreedingEstimator::logProbPGivenGamma(){
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
	double gammaMinusOne = Gamma.getNaturalScaleValue() - 1.0;
	posteriorProbability += gammaMinusOne * sumLogFreq + gammaMinusOne * sumLogOneMinusFreq;
	posteriorProbability += (double) numLoci * (randomGenerator.gammaln(2.0*Gamma.getNaturalScaleValue()) - 2.0*randomGenerator.gammaln(Gamma.getNaturalScaleValue()));

//
//	//new
//	double newP = 0.0;
//	gammaMinusOne = Gamma.getNaturalScaleValue() - 1.0;
//	newP += gammaMinusOne * sumLogFreq + gammaMinusOne * sumLogOneMinusFreq;
//	newP += (double) numLoci * (randomGenerator.gammaln(2.0*Gamma.getNaturalScaleValue()) - 2.0*randomGenerator.gammaln(Gamma.getNaturalScaleValue()));
//	std::cout << "newP " << newP << std::endl;
//
//	//old
//	double old = 0.0;
//	old += (Gamma.getNaturalScaleValue()-1.0) * sumLogFreq + (Gamma.getNaturalScaleValue() - 1.0) * sumLogOneMinusFreq;
//	old += (double) numLoci * (randomGenerator.gammaln(Gamma.getNaturalScaleValue() + Gamma.getNaturalScaleValue()) - randomGenerator.gammaln(Gamma.getNaturalScaleValue()) - randomGenerator.gammaln(Gamma.getNaturalScaleValue()));
//	std::cout << "old " << old << std::endl;



//	if(posteriorProbability > 0)
//		throw "logProbPGivenGamma is larger than 1!";
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
	double trueGamma = 0.2;
	double trueLogGamma = log(trueGamma);
	Gamma.update(trueLogGamma, trueGamma);

	//make grid
	int numProbs = 50;
	float step = 1 / (float) numProbs;
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
			logLikelihood += logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], thisF);
		}

		//add P(p|alpha,beta)
		logLikelihood += logProbPGivenGamma();
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
	double trueGamma = 0.2;
	double trueLogGamma = log(trueGamma);
	Gamma.update(trueLogGamma, trueGamma);

	F.update(0.2, true);


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
		std::cout << "passed params, p[l] "<< p[l] << "newPVlue " << newPValue << " F=" << F.F() << " gamma " << Gamma.getNaturalScaleValue() << std::endl;

		for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
			uint8_t* data = likelihoods.curData();
			if(l == index)
				logLikelihood += logLikelihoodAllInds(data, likelihoods.curSampleSize(), newPValue, F.F());
			if(logLikelihood > 0)
				throw "likelihood is larger than 1!";
		}

		//add P(p|alpha,beta)
		logLikelihood += logProbPGivenGamma();
		outP << newPValue << "\t" << logLikelihood << "\n";
	}
	logfile->done();
}

void TInbreedingEstimator::writeLikelihoodForDebuggingGamma(TParameters & params){
	//open output file
	std::string tracefile = outname + "_logLikelihoodAlpha.txt.gz";
	logfile->list("Will write likelihoods for grid of alpha to '" + tracefile + "'.");
	gz::ogzstream outP(tracefile.c_str());
	if(!outP)
		throw "Failed to open file '" + tracefile + "' for writing!";

	//initialize
	int numProbs = 100;
	float step = 2.0 / (float) numProbs;
	std::cout << "first three true alleleFreq: " << p[0] << " " <<  p[1] << " " << p[2] << std::endl;

	F.update(0.2, true);

	double logLikelihood = 0.0;
	long l = 0;
	for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
		uint8_t* data = likelihoods.curData();
		logLikelihood += logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], F.F());
		if(logLikelihood > 0)
			throw "likelihood is larger than 1!";
	}

	double logLikelihoodAllInds = logLikelihood;

	//calculate ll
	for(int i=0; i<numProbs; ++i){
		double newValue = (double) i*step;
		std::cout << "calculating likelihood for alpha and beta =" << newValue << std::endl;
		if(newValue == 0.0)
			//a small number
			newValue = 0.000000001;
		double newValueLog = log(newValue);

		//only when alpha and beta are equal the beta distribution is symetric and the P(p) = P(1-p),
		// ---->  true MAF should give the same result as the true allele freq
		Gamma.update(newValueLog, newValue);

		//add P(p|alpha,beta)
		logLikelihood = logLikelihoodAllInds + logProbPGivenGamma();
		outP << Gamma.getNaturalScaleValue() << "\t" << logLikelihood << "\n";
	}
	logfile->done();

}

//void TInbreedingEstimator::writeLikelihoodForDebuggingBeta(TParameters & params){
//	//open output file
//	std::string tracefile = outname + "_logLikelihoodBeta.txt.gz";
//	logfile->list("Will write likelihoods for grid of beta to '" + tracefile + "'.");
//	gz::ogzstream outP(tracefile.c_str());
//	if(!outP)
//		throw "Failed to open file '" + tracefile + "' for writing!";
//
//	//initialize
//	double betaValue = 0.0;
//	int numProbs = 100;
//	float step = 2.0 / (float) numProbs;
//	double trueAlpha = 0.5;
//	double trueLogAlpha = log(trueAlpha);
//	alpha.update(trueLogAlpha, trueAlpha);
//	std::cout << "first three true alleleFreq: " << p[0] << " " <<  p[1] << " " << p[2] << std::endl;
//
//	F.update(0.2, true);
//
//	//calculate ll
//	double logLikelihood = 0.0;
//	long l = 0;
//	for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
//		uint8_t* data = likelihoods.curData();
//		logLikelihood += logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], F.F());
//		if(logLikelihood > 0)
//			throw "likelihood is larger than 1!";
//	}
//
//	double logLikelihoodAllInds = logLikelihood;
//
//
//	//add alpha and beta
//	for(int i=0; i<numProbs; ++i){
//		double newBetaValue = betaValue + (double) i*step;
//		std::cout << "calculating likelihood for beta=" << newBetaValue << std::endl;
//		if(newBetaValue == 0.0)
//			//a small number
//			newBetaValue = 0.000000001;
//		double newBetaValueLog = log(newBetaValue);
//		beta.update(newBetaValueLog, newBetaValue);
//
//
//		//add P(p|alpha,beta)
//		logLikelihood = logLikelihoodAllInds + logProbPGivenAlphaBeta();
//		outP << beta.getNaturalScaleValue() << "\t" << logLikelihood << "\n";
//	}
//	logfile->done();
//}

void TInbreedingEstimator::oneMCMCIteration(int iterationNum){
	//update params

	numAcceptedF += updateF();

	long l = 0;
	for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
//		if(l == 3){
			uint8_t* data = likelihoods.curData();
			numAcceptedP[l] += updateP(data, l, likelihoods.curSampleSize(), Gamma);
//		}
	}

	//Gamma
	numAcceptedGamma += updateGamma();
}

void TInbreedingEstimator::printAcceptanceRates(int numIterations){
	logfile->conclude("acceptance rate for F is " + toString((double) numAcceptedF / numIterations));
	logfile->conclude("acceptance rate for gamma is " + toString((double) numAcceptedGamma / numIterations));
	logfile->conclude("acceptance rate for first locus is " + toString((double) numAcceptedP[0] / numIterations));
	logfile->conclude("acceptance rate for second locus is " + toString((double) numAcceptedP[1] / numIterations));
	logfile->conclude("acceptance rate for third locus is " + toString((double) numAcceptedP[2] / numIterations));
	logfile->conclude("acceptance rate for fourth locus is " + toString((double) numAcceptedP[3] / numIterations));
	logfile->conclude("acceptance rate for fifth locus is " + toString((double) numAcceptedP[4] / numIterations));
}


void TInbreedingEstimator::resetAcceptanceRates(){
	numAcceptedF = 0;
	numAcceptedGamma = 0;
	for(unsigned int l=0; l<numLoci; ++l)
		numAcceptedP[l] = 0;
}

void TInbreedingEstimator::adjustProposalWidths(){
	F.adjustProposalWidthAfterBurnin(numAcceptedF, burninLength, logfile);
	p.adjustProposalWidthAfterBurnin(numAcceptedP, burninLength);
	Gamma.adjustProposalWidthAfterBurnin(numAcceptedGamma, burninLength);
}

void TInbreedingEstimator::writeParameterEstimatesOfIteration(std::ofstream & out){
//	out << F.F() << "\t" << Gamma.getNaturalScaleValue() << "\t" << Gamma.getLogValue() << "\t"
//			<< p[0] << "\t" <<  p[1] << "\t" <<  p[2] << "\t" <<  p[3] << "\t" <<  p[117] << std::endl;;

	out << F.F() << "\t" << Gamma.getNaturalScaleValue() << "\t" << Gamma.getLogValue();
	for(int l=0; l<100; ++l){
		out << "\t" << p[l];
	}
	out << "\n";

}

void TInbreedingEstimator::runEstimation(TParameters & params){
	logfile->startIndent("Running MCMC to estimate inbreeding coefficient F and the allele frequency distribution:");

	//open MCMC output file
	std::string tracefile = outname + "_inbreedingMCMC.txt";
	logfile->list("Will write MCMC chain to file '" + tracefile + "'.");
	std::ofstream out(tracefile.c_str());
	if(!out)
		throw "Failed to open file '" + tracefile + "' for writing!";

	tracefile = outname + "_inbreedingMCMC_posteriors.txt.gz";
	logfile->list("Will write posterior distribution statistics to file '" + tracefile + "'.");
	gz::ogzstream outP(tracefile.c_str());
	if(!outP)
		throw "Failed to open file '" + tracefile + "' for writing!";

	//write headers
//	out << "F\talpha\talphaLog\tp[0]\tp[1]\tp[2]\tp[3]\tp[117]\n";
	out << "F\talpha\talphaLog";
	for(int l=0; l<100; ++l){
		out << "\tp[" << l << "]";
	}
	out << "\n";
	outP << "param\tmean_posterior\tvar_posterior\tacceptance_rate\tproposalWidth\ttrue_alleleFreq\n";

	//---------------------------
	//run Burnin(s)
	//---------------------------
	for(int b=0; b<numBurnins; ++b){
		std::string reportString = "Running burnin (" + toString(b) + ") of " + toString(burninLength) + " iterations (";
		logfile->listFlush(reportString + "0%) ...");
		int oldProg=0;
		for(long i=0; i<burninLength; ++i){

			oneMCMCIteration(i);
			writeParameterEstimatesOfIteration(out);

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
		p.setSumsToZero();
	}

	//---------------------------
	//run MCMC
	//---------------------------

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
		}

		// print progress
		int prog = (double) i / (double) numIterations * 100.0;
		if(prog > oldProg){
			oldProg = prog;
			logfile->listOverFlush(progressString + " (" + toString(oldProg) + "%)");
		}
	}
	logfile->overList(progressString + " done!   ");

	//results
	printAcceptanceRates(numIterations);
	outP << "modelWithF\t" << F.posteriorProbModelWithF() << "\t-\t-\t-\t-\n";
	for(unsigned long l=0; l<numLoci; ++l){
		outP << "p[" << l << "]\t" << p.getPosteriorMean(l, numIterations) << "\t" << p.getPosteriorVariance(l, numIterations)<< "\t" << numAcceptedP[l] / (double) numIterations << "\t" << p.proposalWidths[l];
		if(trueAlleleFreqProvided)
			outP << "\t" << trueAlleleFreq[l];
		outP << "\n";
	}

	//clean up
	logfile->endIndent();
	out.close();
	outP.close();

	writeLikelihoodForDebuggingGamma(params);
}
