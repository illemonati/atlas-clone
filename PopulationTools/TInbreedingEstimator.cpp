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
	_logLambda = -1.0;
	_expMinusLambda = -1.0;
	_posteriorProbModelWithF = -1;
}

TInbreedingF::TInbreedingF(double F, float & ProbMovingToModelNoF, double & SdProposal, bool InModelWithF, double Lambda){
	_F = F;
	_probMovingToModelNoF = ProbMovingToModelNoF;
	_sdProposal = SdProposal;
	_inModelWithF = InModelWithF;
	_lambda = Lambda;
	_logLambda = log(Lambda);
	_expMinusLambda = log(1-exp(-_lambda));
	_posteriorProbModelWithF = 0;
}

void TInbreedingF::adjustProposalWidthAfterBurnin(int numAcceptedFModelF, int numIterInModelF){
//	if(numAcceptedF == 0)
//		throw "numAccepted = 0";
	double newProposalWidth = _sdProposal;
	newProposalWidth *=  ((double) numAcceptedFModelF / (double) numIterInModelF) * 3.0;

	if(numIterInModelF == 0){
		_lambda *= 10;
		_logLambda = log(_lambda);
	}

	else if(newProposalWidth / _sdProposal < 0.1)
		newProposalWidth = 0.1 * _sdProposal;
	else if(_sdProposal / newProposalWidth < 0.1)
		newProposalWidth = 10 * _sdProposal;

	else if(newProposalWidth > 1.0)
		newProposalWidth = 1.0;

	else if(newProposalWidth < 0.00001 || numAcceptedFModelF == 0){
		newProposalWidth = 0.00001;
		std::cout << "numAcceptedFModelF is zero!";
	}


	_sdProposal = newProposalWidth;

//	if(!(_sdProposal > 0))
//		throw "standard deviation of F proposal kernel is not larger than 0! It is " + toString(_sdProposal);

}

double TInbreedingF::proposeNew(TRandomGenerator* randomGenerator){
    double newF = _F + randomGenerator->getNormalRandom(0, _sdProposal);
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

void TInbreedingF::updateAndAccept(const double & value, const bool & inModelWithF){
	_F = value;
	_inModelWithF = inModelWithF;
	_posteriorProbModelWithF += inModelWithF;
}

void TInbreedingF::updateAndReject(bool inModelWithF){
	_posteriorProbModelWithF += inModelWithF;
}

void TInbreedingF::resetPosterior(){
	_posteriorProbModelWithF = 0;
}

double TInbreedingF::logPDFExp(const double & thisF){
	return _logLambda - _lambda * thisF;
}

double TInbreedingF::logPDFExp(){
	return logPDFExp(_F);
}

double TInbreedingF::PDFExp(const double & thisF){
	//truncated at 1 -> need to divide by cumulative prob at 1
	return _lambda * exp(-_lambda* thisF) - _expMinusLambda;
}

double TInbreedingF::PDFExp(){
	return PDFExp(_F);
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

double TInbreedingF::proposalWidth(){
	return _sdProposal;
}

//---------------------------
// allele frequencies p
//---------------------------
TAlleleFreq::TAlleleFreq(){
	numLoci = -1;
	numLociModelP = -1;
	lambda = -1.0;
	logLambda = -1.0;
	expMinusLambda = -1.0;
	minAlleleFreq = -1.0;
	probMovingToModel0 = -1.0;
	probMovingToModelP = -1.0;
	logProbMovingToModel0 = -1.0;
	logProbMovingToModelP = -1.0;
}


TAlleleFreq::TAlleleFreq(std::vector<double> & P, double & initialProposalWidthFactor, const int numSamples, double & ProbMovingToModel0, double & ProbMovingToModelP, double & Lambda){
	alleleFreq = P;
	initialAlleleFreq = P;
	numLoci = alleleFreq.size();
	numLociModelP = 0;
	lambda = Lambda;
	logLambda = log(Lambda);
	expMinusLambda = log(1-exp(-lambda));
	minAlleleFreq = 1.0 / ((double) 2.0 * numSamples) / 100.0;
	probMovingToModel0 = ProbMovingToModel0;
	probMovingToModelP = ProbMovingToModelP;
	logProbMovingToModel0 = log(probMovingToModel0);
	logProbMovingToModelP = log(probMovingToModelP);

	//vectors
	//init proposal widths and check range of p
	for(int l=0; l<numLoci; ++l){
		posteriorProbModelP.push_back(0.0);
		if(alleleFreq[l] == 0)
			proposalWidths.push_back(initialProposalWidthFactor * minAlleleFreq);
		else {
			proposalWidths.push_back(initialProposalWidthFactor * alleleFreq[l]);
		}
		sumIterations.push_back(0.0);
		sumOfSquaresIterations.push_back(0.0);

		//set models
		if(alleleFreq[l] == 0.0 || alleleFreq[l] < 0.1*minAlleleFreq){
			modelP.push_back(false);
		} else {
			modelP.push_back(true);
		}
		numLociModelP = numLociModelP + modelP[l];
	}
}

void TAlleleFreq::setSumsForPosteriorToZero(){
	for(int l=0; l<numLoci; ++l){
		posteriorProbModelP[l] = 0.0;
		sumIterations[l] = 0.0;
		sumOfSquaresIterations[l] = 0.0;
	}
}

void TAlleleFreq::setToValue(double fixedValue){
	for(int l=0; l<numLoci; ++l){
		alleleFreq[l] = fixedValue;
	}
}

void TAlleleFreq::adjustProposalWidthAfterBurnin(std::vector<int> & numAcceptedP, std::vector<int> & numUpdates){
	//adjust proposal with for p
	for(long l=0; l<numLoci; ++l){
//		if(modelP[l]){
			if(numAcceptedP[l] == 0){
				std::cout << "!!!!!!!!! acceptance rate is zero for locus " << l << std::endl;
	//			double minFreq = 1.0 / (double) numSamples;
	//			alleleFreq[l] = minFreq;
				proposalWidths[l] *= 0.1;
			} else {
				double newProposalWidth = proposalWidths[l];
				newProposalWidth *=  (double) numAcceptedP[l] / (double) numUpdates[l] * 3.0;

				if(newProposalWidth / proposalWidths[l] < 0.1)
					newProposalWidth = 0.1 * proposalWidths[l];
				else if(proposalWidths[l] / newProposalWidth < 0.1)
					newProposalWidth = 10.0 * proposalWidths[l];

				if(newProposalWidth > 0.5)
					newProposalWidth = 0.5;

				proposalWidths[l] = newProposalWidth;

				if(proposalWidths[l] <= 0)
					throw "Proposal width for allele frequency at locus " + toString(l) + " is not larger than 0!";
			}

//		}
	}
//	numUpdates = 0;
}

void TAlleleFreq::update(long & index, const double & value, const bool ModelP){
	//update numLociModelP
	if(modelP[index] == false && ModelP == true){
		++numLociModelP;
	} else if(modelP[index] == true && ModelP == false){
		--numLociModelP;
	}
	//update current model
	modelP[index] = ModelP;
//	std::cout << "value " << value << std::endl;

	//update posteriors
	posteriorProbModelP[index] += ModelP;
	if(ModelP){
		sumIterations[index] += value;
		sumOfSquaresIterations[index] += value*value;
	}

	//freq value
	if(value < 0 || value > 1)
		throw "updating allele freq at locus " + toString(index) + " to value out of range!";
//	if(value < minAlleleFreq)
//		alleleFreq[index] = minAlleleFreq;
//	else if(value > 1.0 - minAlleleFreq)
//		alleleFreq[index] = 1.0 - minAlleleFreq;
//	else
	alleleFreq[index] = value;

	if(modelP[index] == true && alleleFreq[index] == 0.0)
		throw "allele freq at locus " + toString(index) + " is " + toString(alleleFreq[index]) + " but locus is in modelP=" + toString(modelP[index]);
//	if(modelP[index] == false && alleleFreq[index] > 0.0)
//		throw "allele freq at locus " + toString(index) + " is " + toString(alleleFreq[index]) + " but locus is in modelP=" + toString(modelP[index]);

}

double TAlleleFreq::proposeNew(long & locusNum, TRandomGenerator* randomGenerator){
//	newAlleleFreq[G][l] = fabs(alleleFreq[G][l] + (randomGenerator.getRand()-0.5) * proposalWidthFrequencies[G][l]);

	double newP = (alleleFreq[locusNum] + randomGenerator->getRand() * proposalWidths[locusNum] - proposalWidths[locusNum] / 2.0);
	//is proposed p outside range [0,1]
	if(newP < 0.0){
//		std::cout << "mirroring p at 0: " << newP << " to " << -newP << std::endl;
		newP = - newP;
	} else if(newP > 1.0){
//		std::cout << "mirroring p at 1: " << newP << " to " << 2 - newP << std::endl;
		newP = 2.0 - newP;
	}
	if(newP == 0.0){
		throw "proposing allele frequency zero at locus " + toString(locusNum) + "!";
	}
	return newP;
}

double TAlleleFreq::getPosteriorMean(unsigned long & index, int numUpdates){
	return sumIterations[index] / (double) numUpdates;
}

double TAlleleFreq::getPosteriorVariance(unsigned long & index, int numUpdates){
	return sumOfSquaresIterations[index] / (double) numUpdates - getPosteriorMean(index, numUpdates)*getPosteriorMean(index, numUpdates);
}

double TAlleleFreq::getProposalWidth(const unsigned long & index){
	return proposalWidths[index];
}

long TAlleleFreq::getNumLociInModelP(){
	return numLociModelP;
}

long TAlleleFreq::getNumLociInModel0(){
	return numLoci - numLociModelP;
}

double TAlleleFreq::logPDFExp(const double & thisP){
	//truncated at 1 -> need to divide by cumulative prob at 1
	return logLambda - lambda * thisP - expMinusLambda;
}

double TAlleleFreq::logPDFExp(const long & thisLocus){
	return logPDFExp(alleleFreq[thisLocus]);
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

double TGamma::getProposalWidth(){
	return proposalWidth;
}

//---------------------------
// pi
//---------------------------

TPi::TPi(){
	proposalWidth = -1.0;
	_pi = -1.0;
	_minPi = -1.0;
	_logPi = -1.0;
	_logOneMinusPi = -1.0;
}

TPi::TPi(double & ProposalWidth, double & initialValue){
	proposalWidth = ProposalWidth;
	_pi = initialValue;
	_logPi = log(_pi);
	_logOneMinusPi = log(1.0 - _pi);
	_minPi = 0.000000000000001;
	if(_pi == 0.0)
		_pi = _minPi;
	else if(_pi >= 1)
		_pi = 1 - _minPi;
}

void TPi::update(const double & newNaturalScaleValue){
	if(newNaturalScaleValue < 0)
		throw "updating pi to negative value!";
	_pi = newNaturalScaleValue;
	_logPi = log(_pi);
	_logOneMinusPi = log(1.0 - _pi);
}

void TPi::adjustProposalWidthAfterBurnin(int numAccepted, int numUpdates){
	double newProposalWidth = proposalWidth;
	newProposalWidth *= (double) numAccepted / (double) numUpdates * 3.0;

	if(newProposalWidth / proposalWidth < 0.1)
		newProposalWidth = 0.1 * proposalWidth;
	else if(proposalWidth / newProposalWidth < 0.1)
		newProposalWidth = 10 * proposalWidth;
	if(newProposalWidth > 0.5)
		newProposalWidth = 0.5;

	proposalWidth = newProposalWidth;

	if(proposalWidth <= 0)
		throw "Proposal width for pi is not larger than 0!";
}

double TPi::getPi(){
	return _pi;
}

double TPi::getProposalWidth(){
	return proposalWidth;
}

double TPi::proposeNew(TRandomGenerator* randomGenerator){
	double newPi = _pi + randomGenerator->getRand() * proposalWidth - proposalWidth / 2.0;
	if(newPi < 0.0){
		newPi = -newPi;
	} else if(newPi > 1.0){
		newPi = 2 - newPi;
	}
	return newPi;
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

	//initialize random generator
	//TODO: do the random generator initialization in the task switcher?
	logfile->listFlush("Initializing random generator ...");
	if(Parameters.parameterExists("fixedSeed")){
		randomGenerator = new TRandomGenerator(Parameters.getParameterLong("fixedSeed"), true);
	} else if(Parameters.parameterExists("addToSeed")){
		randomGenerator = new TRandomGenerator(Parameters.getParameterLong("addToSeed"), false);
	} else randomGenerator = new TRandomGenerator();
	logfile->write(" done with seed " + toString(randomGenerator->usedSeed) + "!");

	//initialize parameters
	numAcceptedP.reserve(numLoci);
	numAcceptedPModelP.reserve(numLoci);
	numIterInModelP.reserve(numLoci);
	for(unsigned int i=0; i<numLoci; ++i){
		numAcceptedP.push_back(0);
		numAcceptedPModelP.push_back(0);
		numIterInModelP.push_back(0);
	}

	initParams(randomGenerator, Parameters);
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

	std::cout << "p.numLoci "<< p.numLoci << std::endl;

	for(unsigned int l = 0; l < numLoci; l++){
		if(p.modelP[l]){
			mean += p[l];
			sumXSquare += p[l] * p[l];
		}
	}
	mean /= (double) p.getNumLociInModelP();
	sumXSquare /= (double) p.getNumLociInModelP();
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

	logfile->list("Initialized gamma to " + toString(Gamma.getNaturalScaleValue()));
}

void TInbreedingEstimator::initF(TParameters & parameters){
	sdF = parameters.getParameterDoubleWithDefault("sdF", 0.2);
	logfile->list("Standard deviation of proposal kernel for F is set to " + toString(sdF));

	float probMovingToModelNoF = parameters.getParameterDoubleWithDefault("probMovingToModelNoF", 0.1);
	logfile->list("Will propose move to model without F with probability " + toString(probMovingToModelNoF));

	double lambdaF = parameters.getParameterDoubleWithDefault("lambdaF", 100.0);
	logfile->list("Lambda of exponential distribution used for the proposal of new F after move to model with F is set to " + toString(lambdaF));

	bool startInModelWithF = true;
	if(parameters.parameterExists("startInZeroModel"))
		startInModelWithF = false;
	if(startInModelWithF){
		F = TInbreedingF(randomGenerator->getExponentialRandom(lambdaF), probMovingToModelNoF, sdF, startInModelWithF, lambdaF);
		F.updateAndAccept(randomGenerator->getExponentialRandomTruncated(lambdaF, 0, 1), true);
//		F.updateAndAccept(0.2, true);
		logfile->list("initialized F to " + toString(F.F()) + " in model " + toString(F.inModelWithF()));
	} else {
		F = TInbreedingF(0, probMovingToModelNoF, sdF, startInModelWithF, lambdaF);
		logfile->list("initialized F to " + toString(F.F()) + " in model " + toString(F.inModelWithF()));
	}
}

void TInbreedingEstimator::initAlleleFreq(TParameters & parameters){
	double widthProposalKernelP = parameters.getParameterDoubleWithDefault("widthProposalKernelPFactor", 2.0);
	logfile->list("Will use a proposal kernel of width " + toString(widthProposalKernelP) + " for updates of p within ModelP");

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
	double probMovingToModel0 = parameters.getParameterDoubleWithDefault("probMovingToModelP0", 0.1);
	logfile->list("Will propose move to model with p=0 with probability " + toString(probMovingToModel0));

	double probMovingToModelP = parameters.getParameterDoubleWithDefault("probMovingToModelP", 0.1);
	logfile->list("Will propose move to model with p=0 with probability " + toString(probMovingToModelP));

	double lambdaP = parameters.getParameterDoubleWithDefault("lambdaP", 100.0);
	logfile->list("Lambda of exponential distribution used for the proposal of new p after move to modelP is set to " + toString(lambdaP));

	p = TAlleleFreq(tmp2, widthProposalKernelP, likelihoods.getNumIndividuals(), probMovingToModel0, probMovingToModelP, lambdaP);
//	p.setToValue(0.2);

	if(p.numLoci != numLoci)
		throw "Did not receive one allele frequency per locus! Number of loci=" + toString(numLoci) + " and number of allele frequencies " + toString(p.numLoci);
	logfile->list("initialized allele frequencies for " + toString(p.numLoci) + " loci");

	std::cout << "locus 0 is initialized to " << p[0] << std::endl;
}

void TInbreedingEstimator::initParams(TRandomGenerator* randomGenerator, TParameters & parameters){
	//F
	initF(parameters);

	//p
	initAlleleFreq(parameters);

	//gamma
	double widthProposalKernelGamma = parameters.getParameterDoubleWithDefault("widthProposalKernelGamma", 0.35);
	Gamma = TGamma(widthProposalKernelGamma);

	initializeGamma();
//	Gamma.update(log(0.5), 0.5);

	logfile->list("Will use a proposal kernel of width " + toString(widthProposalKernelGamma) + " for updates of log(alpha)");
	logfile->list("Will use a proposal kernel of width " + toString(widthProposalKernelGamma) + " for updates of log(beta)");

	//pi
	double widthProposalKernelPi = parameters.getParameterDoubleWithDefault("widthProposalKernelPi", 1.0 / (double) p.numLoci);
	double initialValue = (double) p.getNumLociInModelP() / (double) p.numLoci;
//	double initialValue = 0.8;
	pi = TPi(widthProposalKernelPi, initialValue);
	logfile->list("initialized pi to " + toString(pi.getPi()));
}

bool TInbreedingEstimator::updateF(){
	if(F.inModelWithF()){
		//try move to model without F
		if(randomGenerator->getRand() < F.probMovingToModelNoF()){
			double logH = F.logPDFExp() - log(F.probMovingToModelNoF());

			long l = 0;
			for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
				logH += logLikelihoodAllInds(likelihoods.curData(), likelihoods.curSampleSize(), p[l], 0)
				- logLikelihoodAllInds(likelihoods.curData(), likelihoods.curSampleSize(), p[l], F.F());
			}

			//accept?
			if(log(randomGenerator->getRand()) < logH){
				//move to model0
				F.updateAndAccept(0, false);
				return true;
			} else {
				//stay in modelF
				++numIterInModelF;
				F.updateAndReject(true);
				return false;
			}

		} else {
			//propose new F
			++numIterInModelF;
			double newF = F.proposeNew(randomGenerator);

			long l = 0;
			double logH = 0.0;
			for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
				logH += logLikelihoodAllInds(likelihoods.curData(), likelihoods.curSampleSize(), p[l], newF)
				- logLikelihoodAllInds(likelihoods.curData(), likelihoods.curSampleSize(), p[l], F.F());
			}

			//accept?
			if(log(randomGenerator->getRand()) < logH){
				F.updateAndAccept(newF, true);
				++numAcceptedFModelF;
				return true;
			} else {
				//reject new value but still stay in modelF
				F.updateAndReject(true);
				return false;
			}
		}
	}

	//try to move to model with F
	else {
		double newF = randomGenerator->getExponentialRandomTruncated(F.lambda(), 0.0, 1.0);
		double logH = log(F.probMovingToModelNoF()) - F.logPDFExp(newF);

		long l = 0;
		for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
			logH += logLikelihoodAllInds(likelihoods.curData(), likelihoods.curSampleSize(), p[l], newF)
			- logLikelihoodAllInds(likelihoods.curData(), likelihoods.curSampleSize(), p[l], 0);
		}

		//accept?
		if(log(randomGenerator->getRand()) < logH){
			//move to modelF
			F.updateAndAccept(newF, true);
			return true;
		} else {
			//stay in model0
			F.updateAndReject(false);
			return false;
		}
	}
}

bool TInbreedingEstimator::updateP(TSampleLikelihoods* data, long & locusNum, int curSampleSize, TGamma & Gamma){
//	std::cout << "locusNum " << locusNum << std::endl;
	if(p.modelP[locusNum]){
		if(randomGenerator->getRand() < p.probMovingToModel0){
//			std::cout << "proposing move from p to zero ";
			//propose model0
			double gammaNat = Gamma.getNaturalScaleValue();
			double logH = 2.0 * randomGenerator->gammaln(gammaNat)
					+ pi.getLogOneMinusPi()
					- p.logProbMovingToModelP
					+ p.logPDFExp(locusNum)
					- (gammaNat-1.0) * log(p[locusNum])
					- (gammaNat-1.0) * log(1.0-p[locusNum])
					- randomGenerator->gammaln(2.0*gammaNat)
					- pi.getLogPi()
					- p.logProbMovingToModel0 //if prob is zero we never get here
					+ logLikelihoodAllInds(data, likelihoods.curSampleSize(), 0.0, F.F())
					- logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[locusNum], F.F());

			//accept?
			double tmp = log(randomGenerator->getRand());
//			if(locusNum == 26){
//				for(int s=0; s<curSampleSize; ++s){
//					int index = 3*s;
//					//calculate and add ratio for each genotype
//					std::cout << (int) data[index] << "," << (int) data[index + 1] << "," << (int) data[index+2] << " " <<std::endl;
//				}
//
//				std::cout << "curP " << p[26] << " proposed " << "model0" << " hastings " << exp(logH) << " random " << exp(tmp) << std::endl;;
//				std::cout << "ll ratio " << logLikelihoodAllInds(data, likelihoods.curSampleSize(), 0.0, F.F())
//							- logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[locusNum], F.F())
//							<< " pdfExp " <<p.logPDFExp(locusNum) << std::endl;
//			}
			if(tmp < logH){
//				if(locusNum == 26)
//					std::cout << " accepted" << std::endl;
				//move to model0
				p.update(locusNum, 0.0, false);
				return true;
			} else {
				//stay in modelP
//				if(locusNum == 26)
//					std::cout << " rejected" << std::endl;
				++numIterInModelP[locusNum];
				p.update(locusNum, p[locusNum], true);
				return false;
			}
		} else {
			//propose new p
//			std::cout << "proposing move from p to p ";

			++numIterInModelP[locusNum];
			double newP = p.proposeNew(locusNum, randomGenerator);
			if(newP < 0)
				throw "proposed negative newP in move from modelP to modelP': " + toString(newP);

			double logH = (Gamma.getNaturalScaleValue() - 1.0) * (log(newP) + log(1.0 - newP) - log(p[locusNum]) - log(1.0 - p[locusNum])) // - log(1.0 - p[locusNum]));													)
						+ logLikelihoodAllInds(data, likelihoods.curSampleSize(), newP, F.F())
						- logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[locusNum], F.F());

			//accept?
			double tmp = log(randomGenerator->getRand());
//			if(locusNum == 2)
//				std::cout << "curP " << p[2] << " proposed " << newP << " hastings " << logH << " random " << tmp << std::flush;
			if(tmp < logH){
//				if(locusNum == 2 && newP < p[2])
//					std::cout << "curP " << p[2] << " proposed " << newP << " hastings " << exp(logH) << " random " << exp(tmp) << std::endl;
				//update p
				p.update(locusNum, newP, true);
				++numAcceptedPModelP[locusNum];
				return true;
			} else {
				//don't update p
				p.update(locusNum, p[locusNum], true);
				return false;
			}
		}

	} else {
		if(randomGenerator->getRand() < p.probMovingToModelP){
			//move to model p
			double newP = randomGenerator->getExponentialRandomTruncated(p.lambda, 0.0, 1.0);
			if(newP < 0)
				throw "proposed negative newP in move from model0 to modelP: " + toString(newP);
			double gammaNat = Gamma.getNaturalScaleValue();
			double logH = - 2.0 * randomGenerator->gammaln(gammaNat)
					- pi.getLogOneMinusPi()
					+ p.logProbMovingToModelP
					- p.logPDFExp(newP)
					+ (gammaNat-1.0) * log(newP)
					+ (gammaNat-1.0) * log(1.0-newP)
					+ randomGenerator->gammaln(2.0*gammaNat)
					+ pi.getLogPi()
					- logLikelihoodAllInds(data, likelihoods.curSampleSize(), 0.0, F.F())
					+ logLikelihoodAllInds(data, likelihoods.curSampleSize(), newP, F.F());
			if(p.probMovingToModel0 == 0.0){
				logH += log(0.00000000000000000001);
			} else {
				logH += p.logProbMovingToModel0;
			}

			//accept?
			double tmp = log(randomGenerator->getRand());
	//		if(locusNum==0){
	//			std::cout << "model zero to P: logH " << logH << " oldP " << p[locusNum] << " newP " << newP << " log random " << tmp <<" proposalWidth " << p.getProposalWidth(locusNum)<< std::endl;
	//		}

			if(tmp < logH){
	//			std::cout << "accepted" << std::endl;
				p.update(locusNum, newP, true);
				return true;
			} else {
	//			std::cout << "rejected" << std::endl;
				p.update(locusNum, p[locusNum], false);
				return false;
			}

		} else
			return false;
	}
}

bool TInbreedingEstimator::updateGamma(){
	// propose new log(value) (uniform proposal kernel)
	double newLogGamma = Gamma.getLogValue() + randomGenerator->getRand() * Gamma.proposalWidth - Gamma.proposalWidth / 2.0;
	double newGamma = exp(newLogGamma);
	if(newGamma < 0)
		throw "proposed new log(alpha) that on natural scale is negative!";

	//calc sum of log(f) and log(1-f)
	double sumLogFreq = 0.0;
	double sumLogOneMinusFreq = 0.0;

	for(unsigned long l=0; l<numLoci; ++l){
		if(p.modelP[l]){
			sumLogFreq += log(p[l]);
			sumLogOneMinusFreq += log(1.0 - p[l]);
//			if(p[l] == 0.0)
//				logfile->warning("allele freq at locus " + toString(l) +  " is zero, " + toString(p[l]) + " but in modelP=" + toString(p.modelP[l]));
		} else{
//			if(p[l] > 0.0)
//				logfile->warning("allele freq at locus " + toString(l) +  " is not zero, " + toString(p[l]) + " but in modelP=" + toString(p.modelP[l]));
		}
	}

	// compute log hastings ratio
	double diffGammas = newGamma - Gamma.getNaturalScaleValue();
	double logH = (double) p.getNumLociInModelP() *
			(randomGenerator->gammaln(2.0*newGamma)
			+ 2.0*randomGenerator->gammaln(Gamma.getNaturalScaleValue())
			- randomGenerator->gammaln(2.0*Gamma.getNaturalScaleValue())
			- 2.0*randomGenerator->gammaln(newGamma)
			)
			+ diffGammas * sumLogFreq
			+ diffGammas * sumLogOneMinusFreq;

	// accept or reject
	double tmp = log(randomGenerator->getRand());
	if(tmp < logH){
		Gamma.update(newLogGamma, newGamma);
		return true;
	}
	else
		return false;
}

bool TInbreedingEstimator::updatePi(){
	// propose new log(value) (uniform proposal kernel)
	double newPi = pi.proposeNew(randomGenerator);
	double numLociInModelP = (double) p.getNumLociInModelP();
	double numLociInModel0 = (double) p.getNumLociInModel0();

	double logH = numLociInModelP * log(newPi) + numLociInModel0 * log(1.0-newPi)
			 - numLociInModelP * log(pi.getPi()) - numLociInModel0 * log(1.0-pi.getPi());

	// accept or reject
	double tmp = log(randomGenerator->getRand());
//	std::cout << " curPi " << pi.getPi() << " newPi " << newPi << " numLociInModelP " << numLociInModelP  << " numLociInModel0 " << numLociInModel0<< " H " << exp(logH) << " tmp " << exp(tmp)<< std::flush;

	if(tmp < logH){
		pi.update(newPi);
//		std::cout << " accepted" << std::endl;
		return true;
	}
	else{
//		std::cout << " rejected" << std::endl;
		return false;
	}
}

double TInbreedingEstimator::logLikelihoodAllInds(TSampleLikelihoods* data, int curSampleSize, double thisP, double thisF){
	if(thisP < 0)
		throw "allele freq is negative!";
	//sum over all individuals of log sum_g P(d|g)P(g|p,F)
	double sumOverInds = 0.0;

	double PGeno[3];
	if(thisP == 0.0){
		PGeno[0] = 1.0;
		PGeno[1] = 0.0;
		PGeno[2] = 0.0;
	} else {
		PGeno[0] = (1.0 - thisF)*(1.0 - thisP)*(1.0 - thisP) + thisF*(1.0 - thisP);
		PGeno[1] = 2.0*thisP*(1.0 - thisP)*(1.0 - thisF);
		PGeno[2] = (1.0 - thisF)*thisP*thisP + thisF*thisP;
	}

	for(int s=0; s<curSampleSize; ++s){
		if(!data[s].isMissing && !data[s].isHaploid){
			//calculate and add ratio for each genotype
			double integrationOverGeno = qualMap.phredToError(data[s].phredLikelihood_0) * PGeno[0];
			integrationOverGeno += qualMap.phredToError(data[s].phredLikelihood_1) * PGeno[1];
			integrationOverGeno += qualMap.phredToError(data[s].phredLikelihood_2) * PGeno[2];

			//check if likelihood of sample is a probability
			if(integrationOverGeno < 0){
				throw "Probability of genotype is negative: " + toString(integrationOverGeno);
			}

			sumOverInds += log(integrationOverGeno);
		}
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
//		if(l == 736){
			if(p.modelP[l]){
				if(p[l] == 0.0)
					throw "allele freq is zero when taking log";
				sumLogFreq += log(p[l]);
				sumLogOneMinusFreq += log(1.0 - p[l]);
			}
//		}
	}

	//add beta density
	double gammaMinusOne = Gamma.getNaturalScaleValue() - 1.0;
	posteriorProbability += gammaMinusOne * sumLogFreq + gammaMinusOne * sumLogOneMinusFreq;
	posteriorProbability += (double) p.getNumLociInModelP() * (randomGenerator->gammaln(2.0*Gamma.getNaturalScaleValue()) - 2.0*randomGenerator->gammaln(Gamma.getNaturalScaleValue()));

//	std::cout << "gamma " << Gamma.getNaturalScaleValue() << " F " << F.F() << " p[26] " << p[26] << " modelP " << p.modelP[26] << " pi " << pi.getPi()<< "posteriorProbability " << posteriorProbability << std::endl;

	return posteriorProbability;
};

double TInbreedingEstimator::getLogLikelihoodCurrentParams(){
	double logLikelihood = 0.0;
	long l = 0;
	for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
//		if(l == 736){
			logLikelihood += logLikelihoodAllInds(likelihoods.curData(), likelihoods.curSampleSize(), p[l], F.F());
//		}
	}

	//add P(p|alpha,beta)
	logLikelihood += logProbPGivenGamma();
//	std::cout << "logLikelihood " << logLikelihood << std::endl;
	return logLikelihood;
};

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
			logLikelihood += logLikelihoodAllInds(likelihoods.curData(), likelihoods.curSampleSize(), p[l], thisF);
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
	long index = params.getParameterInt("locus");

	//open output file
	std::string tracefile = outname + "_logLikelihoodP[" + toString(index) + "].txt.gz";
	logfile->list("Will write likelihood chain for grid of p[" + toString(index) + "] to file '" + tracefile + "'.");
	logfile->list("The true value for p[" + toString(index) + "] is: " + toString(p[index]));
	gz::ogzstream outP(tracefile.c_str());
	if(!outP)
		throw "Failed to open file '" + tracefile + "' for writing!";
	std::cout << "first three true alleleFreq: " << p[0] << " " <<  p[1] << " " << p[2] << std::endl;

	//initialize
	double trueGamma = 0.5;
	double trueLogGamma = log(trueGamma);
	Gamma.update(trueLogGamma, trueGamma);

	F.updateAndAccept(0.2, true);
	pi.update(0.8);


	//make grid
	int numProbs = 1000;
	float step = 1.0 / (float) numProbs;
	double pValue = 0.0;

	//calculate ll
	for(int i=0; i<numProbs; ++i){
		double newPValue = pValue + (double) i*step;
		std::cout << "calculating likelihood for p[" << index << "]=" << newPValue << std::endl;
		if(newPValue == 0)
			p.update(index, newPValue, false);
		else
			p.update(index, newPValue, true);


		double logLikelihood = getLogLikelihoodCurrentParams();
		outP << newPValue << "\t" << logLikelihood << "\n";
	}
	logfile->done();
}

void TInbreedingEstimator::writeLikelihoodForDebuggingGamma(TParameters & params){
	//open output file
	std::string tracefile = outname + "_logLikelihoodGamma.txt.gz";
	logfile->list("Will write likelihoods for grid of gamma to '" + tracefile + "'.");
	gz::ogzstream outP(tracefile.c_str());
	if(!outP)
		throw "Failed to open file '" + tracefile + "' for writing!";

	//initialize
	int numProbs = 100;
	float step = 2.0 / (float) numProbs;
	std::cout << "first three true alleleFreq: " << p[0] << " " <<  p[1] << " " << p[2] << std::endl;

	F.updateAndAccept(0.2, true);

	double logLikelihood = 0.0;
	long l = 0;
	for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
		logLikelihood += logLikelihoodAllInds(likelihoods.curData(), likelihoods.curSampleSize(), p[l], F.F());
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

//void TInbreedingEstimator::writeLikelihoodForDebuggingPi(TParameters & params){
//	//open output file
//	std::string tracefile = outname + "_logLikelihoodPi.txt.gz";
//	logfile->list("Will write likelihoods for grid of pi to '" + tracefile + "'.");
//	gz::ogzstream outP(tracefile.c_str());
//	if(!outP)
//		throw "Failed to open file '" + tracefile + "' for writing!";
//
//	//initialize
//	int numProbs = 100;
//	float step = 2.0 / (float) numProbs;
//	std::cout << "first three true alleleFreq: " << p[0] << " " <<  p[1] << " " << p[2] << std::endl;
//
//	F.updateAndAccept(0.2, true);
//	Gamma.update(log(0.5), 0.5);
//
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
//	for(int i=0; i<numProbs; ++i){
//		double newValue = (double) i*step;
//		pi.update(newValue);
//	}
//
//}

void TInbreedingEstimator::oneMCMCIteration(int iterationNum){
	//update params

	numAcceptedF += updateF();

	long l = 0;
	for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
//		if(l == 736){
			numAcceptedP[l] += updateP(likelihoods.curData(), l, likelihoods.curSampleSize(), Gamma);
//		}
	}

	//Gamma
	numAcceptedGamma += updateGamma();
	numAcceptedPi += updatePi();
}

void TInbreedingEstimator::printAcceptanceRates(int numIterations){
	numIterations = (double) numIterations;
	std::cout << "numAcceptedPModelP.at(0) "  << numAcceptedPModelP.at(0) << std::endl;
	std::cout << "numIterInModelP[0] "  << numIterInModelP[0] << std::endl;

	logfile->conclude("total acceptance rate for F is " + toString((double) numAcceptedF / numIterations) + " and acceptance rate for moves within M_F " + toString((double) numAcceptedFModelF / (double) numIterInModelF));
	logfile->conclude("total acceptance rate for gamma is " + toString((double) numAcceptedGamma / numIterations));
	logfile->conclude("total acceptance rate for pi is " + toString((double) numAcceptedPi / numIterations));
	logfile->conclude("total acceptance rate for first locus is " + toString((double) numAcceptedP.at(0) / numIterations) + " and acceptance rate for moves within M_p is " + toString((double) numAcceptedPModelP.at(0) / numIterInModelP[0]));
	logfile->conclude("total acceptance rate for second locus is " + toString((double) numAcceptedP.at(1) / numIterations) + " and acceptance rate for moves within M_p is " + toString((double) numAcceptedPModelP.at(1) / numIterInModelP[1]));
	logfile->conclude("total acceptance rate for third locus is " + toString((double) numAcceptedP.at(2) / numIterations) + " and acceptance rate for moves within M_p is " + toString((double) numAcceptedPModelP.at(2) / numIterInModelP[2]));
	logfile->conclude("total acceptance rate for fourth locus is " + toString((double) numAcceptedP.at(3) / numIterations) + " and acceptance rate for moves within M_p is " + toString((double) numAcceptedPModelP.at(3) / numIterInModelP[3]));
	logfile->conclude("total acceptance rate for fifth locus is " + toString((double) numAcceptedP.at(4) / numIterations) + " and acceptance rate for moves within M_p is " + toString((double) numAcceptedPModelP.at(4) / numIterInModelP[4]));
}


void TInbreedingEstimator::resetAcceptanceRates(){
	numAcceptedF = 0;
	numAcceptedFModelF = 0;
	numIterInModelF = 0;
	F.resetPosterior();
	numAcceptedGamma = 0;
	for(unsigned int l=0; l<numLoci; ++l){
		numAcceptedP[l] = 0;
		numAcceptedPModelP[l] = 0;
		numIterInModelP[l] = 0;
	}
	p.setSumsForPosteriorToZero();
	numAcceptedPi = 0;
}

void TInbreedingEstimator::adjustProposalWidths(){
	F.adjustProposalWidthAfterBurnin(numAcceptedFModelF, numIterInModelF);
	p.adjustProposalWidthAfterBurnin(numAcceptedPModelP, numIterInModelP);
	Gamma.adjustProposalWidthAfterBurnin(numAcceptedGamma, burninLength);
	pi.adjustProposalWidthAfterBurnin(numAcceptedPi, burninLength);
}

void TInbreedingEstimator::writeParameterEstimatesOfIteration(std::ofstream & out){
//	out << F.F() << "\t" << Gamma.getNaturalScaleValue() << "\t" << Gamma.getLogValue() << "\t"
//			<< p[0] << "\t" <<  p[1] << "\t" <<  p[2] << "\t" <<  p[3] << "\t" <<  p[117] << std::endl;;

	out << F.F() << "\t" << Gamma.getNaturalScaleValue() << "\t" << Gamma.getLogValue() << "\t" << pi.getPi();
	for(int l=0; l<100; ++l){
		out << "\t" << p[l];
	}

	out << "\t" << getLogLikelihoodCurrentParams();
	out << "\t" << p.getNumLociInModelP();
//	out << "\t" << numAcceptedPModelP[0] << "\t" << numIterInModelP[0];
	out << "\n";

/*
	//-------------------------------------
	// DEBUG F
	//-------------------------------------



	//M_F -> M_0
	double logH = F.logPDFExp(0.05) - log(F.probMovingToModelNoF());

	out << "\t" << logH;

	long l = 0;
	for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
		uint8_t* data = likelihoods.curData();
		logH += logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], 0)
		- logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], 0.05);
	}

	out << "\t" << logH;


	//M_0 -> M_F
	double newF = 0.05;
	logH = log(F.probMovingToModelNoF()) - F.logPDFExp(newF);

	out << "\t" << logH;

	l = 0;
	for(likelihoods.begin(); !likelihoods.end(); likelihoods.next(), ++l){
		uint8_t* data = likelihoods.curData();
		logH += logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], newF)
		- logLikelihoodAllInds(data, likelihoods.curSampleSize(), p[l], 0);
	}

	out << "\t" << logH;


	out << "\n";
*/
	//-------------------------------------
	// DEBUG p
	//-------------------------------------

//	//propose model0
//	likelihoods.begin();
//	uint8_t* data = likelihoods.curData();
//	double gammaNat = Gamma.getNaturalScaleValue();
//	double logH = 2.0 * randomGenerator->gammaln(gammaNat)
//			+ log(1.0 - pi.getPi())
//			+ p.logPDFExp(0.2)
//			- (gammaNat-1.0) * log(0.2)
//			- (gammaNat-1.0) * log(1.0-0.2)
//			- randomGenerator->gammaln(2.0*gammaNat)
//			- log(pi.getPi())
//			- log(p.probMovingToModel0)
//			+ logLikelihoodAllInds(data, likelihoods.curSampleSize(), 0.0, F.F())
//			- logLikelihoodAllInds(data, likelihoods.curSampleSize(), 0.2, F.F());
//
//	std::cout << "\n--------------" << std::endl;
//
//	std::cout << "p to 0: " << logH << std::endl;
//
//	double newP = 0.2;
//	if(newP < 0)
//		throw "proposed negative newP in move from model0 to modelP: " + toString(newP);
//	logH = - 2.0 * randomGenerator->gammaln(gammaNat)
//			- log(1.0 - pi.getPi())
//			- p.logPDFExp(newP)
//			+ (gammaNat-1.0) * log(newP)
//			+ (gammaNat-1.0) * log(1.0-newP)
//			+ randomGenerator->gammaln(2.0*gammaNat)
//			+ log(pi.getPi())
//			+ log(p.probMovingToModel0)
//			- logLikelihoodAllInds(data, likelihoods.curSampleSize(), 0.0, F.F())
//			+ logLikelihoodAllInds(data, likelihoods.curSampleSize(), newP, F.F());
//
//	std::cout << "0 to p: " << logH << std::endl;

}

void TInbreedingEstimator::writePosteriors(int i){
	std::string tracefile = outname + "_" + toString(i) + "_inbreedingMCMC_posteriors.txt.gz";
	logfile->list("Will write posterior distribution statistics to file '" + tracefile + "'.");
	gz::ogzstream outP(tracefile.c_str());
	if(!outP)
		throw "Failed to open file '" + tracefile + "' for writing!";

	outP << "param\tprop_in_nonZeroModel\tmean_posterior\tvar_posterior\tacceptance_rate_modelNonZero\tproposalWidth\ttrue_alleleFreq\n";
	outP << "F\t" << (double) F.posteriorProbModelWithF() / (double) numIterations  << "\tNA\tNA\t" << numAcceptedFModelF<< "\t" << F.proposalWidth() << "\tNA\n";
	for(unsigned long l=0; l<numLoci; ++l){
		outP << "p[" << l << "]\t" << p.posteriorProbModelP[l] << "\t"<< p.getPosteriorMean(l, numIterations) << "\t" << p.getPosteriorVariance(l, numIterations)<< "\t" << numAcceptedPModelP[l] / (double) numIterInModelP[l] << "\t" << p.getProposalWidth(l);
		if(trueAlleleFreqProvided)
			outP << "\t" << trueAlleleFreq[l];
		else
//			outP << "\tNA";
			outP << "\t" << p.initialAlleleFreq[l];
		outP << "\n";
	}
	outP.close();

}

void TInbreedingEstimator::runEstimation(TParameters & params){
	logfile->startIndent("Running MCMC to estimate inbreeding coefficient F and the allele frequency distribution:");

	//open MCMC output file
	std::string tracefile = outname + "_inbreedingMCMC.txt";
	logfile->list("Will write MCMC chain to file '" + tracefile + "'.");
	std::ofstream out(tracefile.c_str());
	if(!out)
		throw "Failed to open file '" + tracefile + "' for writing!";

	//write header of trace file
	out << "F\tgamma\tgammaLog\tpi";
	for(int l=0; l<100; ++l){
		out << "\tp[" << l << "]";
	}

//	//DEBUG
//	out << "\tq(FtoZero)/P(F')\tH_M0_to_MF\tq(P(F)/q(FtoZero)\tH_MF_to_M0";


	out << "\tlogLikelihood\tp.getNumLociInModelP()\n";
	writeParameterEstimatesOfIteration(out);


	//---------------------------
	//run Burnin(s)
	//---------------------------
	for(int b=0; b<numBurnins; ++b){
		std::string reportString = "Running burnin (" + toString(b) + ") of " + toString(burninLength) + " iterations (";
		logfile->listFlush(reportString + "0%) ...");
		int oldProg=0;
		for(long i=0; i<burninLength; ++i){

			oneMCMCIteration(i);
			if(params.parameterExists("writeBurninToFile"))
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

		//report new proposal widths
		logfile->conclude("New proposal width for F is " + toString(F.proposalWidth()));
		for(int l=0; l<5; ++l)
			logfile->conclude("New proposal width for allele freq at locus one " + toString(p.getProposalWidth(l)));
		logfile->conclude("New proposal width for Gamma is " + toString(Gamma.getProposalWidth()));
		logfile->conclude("New proposal width for pi is " + toString(pi.getProposalWidth()));

		resetAcceptanceRates();
		p.setSumsForPosteriorToZero();
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

		if(i % 1000 == 0){
			writePosteriors(i);
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
	writePosteriors(numIterations);

	//clean up
	logfile->endIndent();
	out.close();

//	writeLikelihoodForDebuggingGamma(params);
}
