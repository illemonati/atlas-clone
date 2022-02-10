/*
 * TThetaEstimator.cpp
 *
 *  Created on: Sep 25, 2017
 *      Author: phaentu
 */

#include "TThetaEstimator.h"
#include "GenotypeTypes.h"

namespace GenotypeLikelihoods{

//---------------------------------------------------------------
//TThetaEstimator_base
//---------------------------------------------------------------
TThetaEstimator_base::TThetaEstimator_base(coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator){
	logfile = Logfile;
	randomGenerator = RandomGenerator;
	useTmpFile = false;
	minSitesWithData = 0;
	data = NULL;
	dataInitialized = false;
	extraVerbose = false;

	initialTheta = 0.01;
	initThetaNumSearchIterations = 10;
	initThetaSearchFactor = 100;
};

TThetaEstimator_base::TThetaEstimator_base(coretools::TParameters & params, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator){
	logfile = Logfile;
	randomGenerator = RandomGenerator;
	data = NULL;
	dataInitialized = false;

	//data
	useTmpFile = params.parameterExists("useTmpFile");
	initDataStorage();

	//minimum window size
	minSitesWithData = params.getParameterWithDefault<int>("minSitesWithData", 1000);

	//extra verbosity
	extraVerbose = params.parameterExists("extraVerbose");
};

TThetaEstimator_base::TThetaEstimator_base(const TThetaEstimator_base & other){
	logfile = other.logfile;
	randomGenerator = other.randomGenerator;
	data = NULL;
	dataInitialized = false;

	//data
	useTmpFile = other.useTmpFile;
	initDataStorage();

	//minimum window size
	minSitesWithData = other.minSitesWithData;

	//extra verbosity
	extraVerbose = other.extraVerbose;
};

void TThetaEstimator_base::initDataStorage(){
	if(dataInitialized){
		delete data;
	}

	//file or vector?
	if(useTmpFile){
		tmpFileName = "temporaryDataForThetaEstimation_" + randomGenerator->getRandomAlphaNumericString(10) + ".tmp.gz";
		logfile->list("Will write temporary data to file '" + tmpFileName + "'.");
		data = new TThetaEstimatorDataFile(tmpFileName);
	} else {
		data = new TThetaEstimatorDataVector();
	}
	dataInitialized = true;
};

void TThetaEstimator_base::readParametersRegardingInitialSearch(coretools::TParameters & params){
	logfile->startIndent("Parameters of the initial theta search:");
	initialTheta = params.getParameterWithDefault("initTheta", 0.01);
	logfile->list("Will start with an initial theta of " + toString(initialTheta) + ".");
	initThetaNumSearchIterations = params.getParameterWithDefault("initThetaNumSearchIterations", 10);
	if(initThetaNumSearchIterations > 0){
		logfile->list("Will run " + toString(initThetaNumSearchIterations) + " iterations of a crude search for an initial theta.");
		initThetaSearchFactor = params.getParameterWithDefault("initThetaSearchFactor", 100);
		logfile->list("The initial search factor will be " + toString(initThetaSearchFactor) + ".");
	} else {
		initThetaSearchFactor = 0;
	}
	logfile->endIndent();
};

void TThetaEstimator_base::fillPGenotype(TGenotypeProbabilities & pGeno, double expTheta, const TBaseProbabilities & baseFrequencies){
	using namespace genometools;
	//assumes that base frequencies are set!
	for(Base b = Base::min; b < Base::max; ++b){
		//homozygous genotypes
		Genotype hom = genotype(b, b);
		pGeno[hom] = baseFrequencies[b] * (expTheta + baseFrequencies[b].get() * (1.0 - expTheta));
		//heterozygous genotypes
		for(Base c = next(b); c < Base::max; ++c){
			Genotype het = genotype(b, c);
			pGeno[het] = 2.0 * baseFrequencies[b].get() * baseFrequencies[c].get() *  (1.0 - expTheta);
		}
	}
};

void TThetaEstimator_base::fillPGenotype(GenotypeLikelihoods::TGenotypeProbabilities & pGeno, const Theta & thisTheta){
	fillPGenotype(pGeno, thisTheta.expTheta, thisTheta.baseFreq);
};

void TThetaEstimator_base::findGoodStartingTheta(TThetaEstimatorData* thisData, Theta & thisTheta, std::string tag){
	logfile->listFlush("Estimating initial parameters" + tag + " ...");

	//set base frequencies to initial base frequencies
	thisTheta.baseFreq = thisData->baseFrequencies();

	//variables
	double initTheta = initialTheta;
	double oldTheta = initTheta;
	double expTheta = exp(-initTheta);

	//calc initial LL
	fillPGenotype(pGenotype, expTheta, thisTheta.baseFreq);
	thisTheta.LL = thisData->calcLogLikelihood(pGenotype);

	//run iterations
	double oldLL = thisTheta.LL;
	double factor = initThetaSearchFactor;
	int numUpdates;
	for(int i=0; i<initThetaNumSearchIterations; ++i){
		//first test increase in theta
		numUpdates = -1;
		do{
			++numUpdates;
			oldLL = thisTheta.LL;
			oldTheta = initTheta;
			initTheta *= factor;
			expTheta = exp(-initTheta);
			fillPGenotype(pGenotype, expTheta, thisTheta.baseFreq);
			thisTheta.LL = thisData->calcLogLikelihood(pGenotype);
		} while(oldLL < thisTheta.LL);
		if(numUpdates == 0){
			//then test decrease in theta
			initTheta = oldTheta;
			thisTheta.LL = oldLL;
			//maybe smaller?
			do{
				oldLL = thisTheta.LL;
				oldTheta = initTheta;
				initTheta /= factor;
				expTheta = exp(-initTheta);
				fillPGenotype(pGenotype, expTheta, thisTheta.baseFreq);
				thisTheta.LL = thisData->calcLogLikelihood(pGenotype);
			} while(oldLL < thisTheta.LL);
		}
		factor = sqrt(factor);
		initTheta = oldTheta;
		thisTheta.LL = oldLL;
	}
	//return previous
	thisTheta.setTheta(oldTheta);
	thisTheta.LL = oldLL;

	//check if values make sense. If theta < 1/(10*windowsize), set it to 1/(10*windowsize)
	if(thisTheta.theta < 0.1/(double) thisData->size()){
		thisTheta.setTheta(0.1/(double) thisData->size());
	} else if(thisTheta.theta > 1.0){
		thisTheta.setTheta(1.0);
	}

	//conclude
	logfile->done();
	logfile->conclude("Initial base frequencies: " + thisTheta.getBaseFrequencyString());
	logfile->conclude("Initial theta = ", thisTheta.theta);
}

//-------------------------------------------------------
//TThetaEstimator
//-------------------------------------------------------
TThetaEstimator::TThetaEstimator(TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator):TThetaEstimator_base(params, Logfile, RandomGenerator){
	estimationSuccessful = false;

	//parse
	logfile->startIndent("Parameters of EM algorithm to infer theta:");
	numIterations = params.getParameterWithDefault<int>("iterations", 100);
	logfile->list("Will run up to " + toString(numIterations) + " iterations.");
	numThetaOnlyUpdates = params.getParameterWithDefault<int>("iterationsThetaOnly", 10);
	logfile->list("In each iteration, theta will be updated " + toString(numThetaOnlyUpdates) + " times.");

	maxEpsilon = params.getParameterWithDefault("maxEps", 0.000001);
	logfile->list("Will run EM until deltaLL < " + toString(maxEpsilon) + ".");
	NewtonRaphsonNumIterations = params.getParameterWithDefault<int>("NRiterations", 10);
	logfile->list("Will run Newton-Raphson algorithm up to " + toString(NewtonRaphsonNumIterations) + " times.");
	NewtonRaphsonMaxF = params.getParameterWithDefault("maxF", 0.00001);
	logfile->list("Will run Newton-Raphson algorithm until max(F) < " + toString(NewtonRaphsonMaxF) + ".");
	logfile->endIndent();

	//params regarding initial search
	readParametersRegardingInitialSearch(params);
}

TThetaEstimator::TThetaEstimator(TLog* Logfile, TRandomGenerator* RandomGenerator):TThetaEstimator_base(Logfile, RandomGenerator){
	//set EM parameters to default
	numIterations = -1;
	numThetaOnlyUpdates = -1;
	maxEpsilon = 0.0;
	NewtonRaphsonNumIterations = -1;
	NewtonRaphsonMaxF = 0.0;
	initialTheta = 0.0;
	initThetaSearchFactor = -1;
	initThetaNumSearchIterations = -1;
	estimationSuccessful = false;
};

TThetaEstimator::TThetaEstimator(const TThetaEstimator & other):TThetaEstimator_base(other){
	//set EM parameters to default
	numIterations = other.numIterations;
	numThetaOnlyUpdates = -other.numThetaOnlyUpdates;
	maxEpsilon = other.maxEpsilon;
	NewtonRaphsonNumIterations = other.NewtonRaphsonNumIterations;
	NewtonRaphsonMaxF = other.NewtonRaphsonMaxF;
	initialTheta = other.initialTheta;
	initThetaSearchFactor = other.initThetaSearchFactor;
	initThetaNumSearchIterations = other.initThetaNumSearchIterations;
	estimationSuccessful = other.estimationSuccessful;
};

void TThetaEstimator::clear(){
	data->clear();
	estimationSuccessful = false;
};

void TThetaEstimator::add(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods){
	data->add(site, genotypeLikelihoods);
};

void TThetaEstimator::add(const TWindow & window, const TGenotypeLikelihoodCalculator & glCalculator){
	TGenotypeLikelihoods genoLik;
	for(std::vector<TSite>::const_iterator it = window.cbegin(); it != window.cend(); ++it){
		glCalculator.calculateGenotypeLikelihoods(*it, genoLik);
		add(*it, genoLik);
	}
};

double TThetaEstimator::_calcFisherInfo(const TGenotypeProbabilities & _pGenotype, const TGenotypeData deriv_pGenotype){
//sum Ri over all sites
	double FisherInfo = 0.0;
	data->begin();
	do{
		//calc Ri
		double Ri_a = data->curGenotypeLikelihoods().weightedSum(deriv_pGenotype);
		double Ri_b = data->curGenotypeLikelihoods().weightedSum(_pGenotype);
		double Ri = Ri_a / Ri_b;

		//add to Fisher Info
		FisherInfo += Ri * (Ri + 1.0);
	} while(data->next());

	return(FisherInfo);
};

bool TThetaEstimator::_NRAllParams(){
	using namespace genometools;
	//calculate substitution probabilities
	fillPGenotype(pGenotype, theta);

	//Calculate all genotype probabilities for all sites
	data->fillP_G(P_G, pGenotype);

	//prepare storage
	arma::mat Jacobian(6,6);
	arma::vec F(6);
	arma::mat JxF(6,6);

	double rho = theta.expTheta / (1.0 - theta.expTheta);
	double mu = data->sizeWithData();

	TBaseProbabilities& baseFreq = theta.baseFreq; //for cleaner code
	for(int n=0; n<NewtonRaphsonNumIterations; ++n){
		//i) calculate F (Note: index is zero based!)
		F(4) = data->sizeWithData();
		F(5) = 0.0;
		for(Base k = Base::min; k < Base::max; ++k){
			const auto hom = genotype(k, k);
			double tmpSum = 0.0;
			for(Base l = Base::min; l < Base::max; ++l){
				if(l != k){
					const auto het = genotype(k, l);
					tmpSum += P_G[het];
				}
			}
			F(index(k)) = P_G[hom] * (1.0 + baseFreq[k].get() / (rho + baseFreq[k].get())) + tmpSum - mu * baseFreq[k].get();
			F(4) -= P_G[hom] * (rho + 1.0 ) / (rho + baseFreq[k].get());
			F(5) += baseFreq[k].get();
		}
		F(5) = F(5) - 1.0;

		//ii) fill Jacobian (Note: index is zero based!)
		Jacobian.zeros();
		TBaseData tmp;
		for(Base k = Base::min; k < Base::max; ++k){
			const auto hom = genotype(Base(k), Base(k));
			tmp[k] = P_G[hom] / ((baseFreq[k].get() + rho)*(baseFreq[k].get() + rho));
		}

		for(Base k = Base::min; k < Base::max; ++k){
			const auto i  = index(k);
			Jacobian(i,i) = tmp[k] * rho - mu;
			Jacobian(i,4) = - tmp[k];
			Jacobian(5,i) = 1.0;
			Jacobian(4,i) = tmp[k] * (rho + 1.0);
			Jacobian(i,5) = - baseFreq[k].get();
			Jacobian(4,4) += tmp[k] * (1.0 - baseFreq[k].get());
		}

		//iii) now estimate new parameters
		mu = data->sizeWithData();

		if(solve(JxF, Jacobian, F)){
			for(Base k = Base::min; k < Base::max; ++k){
				baseFreq[k] = baseFreq[k].get() - JxF(index(k));
			}
			rho -= JxF(4);
			mu -= JxF(5);

			//check if we break
			double maxF = 0.0;
			for(int i=0; i<6; ++i){
				if(F(i) > maxF) maxF = F(i);
			}

			if(maxF < NewtonRaphsonMaxF){
				theta.setTheta(-log(rho / (1.0 + rho)));
				return true;
			}
		} else {
			return false;
		}
	}
	return true;
};

void TThetaEstimator::_NROnlyTheta(){
	using namespace genometools;
	//calculate	substitution probabilities
	fillPGenotype(pGenotype, theta);

	//Calculate all genotype probabilities for all sites
	data->fillP_G(P_G, pGenotype);

	double rho = theta.expTheta / (1.0 - theta.expTheta);

	for(int n=0; n<NewtonRaphsonNumIterations; ++n){

		//i) calculate F() (Note: index is zero based!)
		//ii) fill Jacobian (Note: index is zero based!)
		double F = data->sizeWithData();
		double Jacobian = 0.0;
		for(Base k = Base::min; k < Base::max; ++k){
			const auto hom = genotype(k, k);

			const double tmp = (theta.baseFreq[k].get() + rho);
			F -= P_G[hom] * (rho + 1.0 ) / tmp;

			const double tmpSum = P_G[hom] / (tmp*tmp);
			Jacobian += tmpSum * (1.0 - theta.baseFreq[k].get());
		}

		//iii) now estimate new parameters
		rho = rho - F / Jacobian;

		//check if we break
		if(F < NewtonRaphsonMaxF){
			theta.setTheta(-log(rho / (1.0 + rho)));
			break;
		}
	}
};

void TThetaEstimator::_runEMForTheta(){
	//increase initialTheta if we fail to calculate inverse of Jacobian.
	// this may be the case if initialTheta is smaller than true theta and likelihood is very flat
	int failedAttempts = 0;
	double startingTheta = initialTheta;
	theta.LL = -9e100;
	while(!_NRAllParams()){
		++failedAttempts;
		//solve did not work -> start with higher theta!
		startingTheta *= 2.0;
		theta.setTheta(startingTheta);
		if(startingTheta > 1.0)
			throw "Failed to estimate Theta, issues calculating inverse of Jacobian!";
	}

	//start EM loop
	for(int iter = 0; iter < numIterations; ++iter){
		//update only theta: most difficult parameter and it is much faster to update only this one alone.
		int i=0;
		double oldTheta = 0.0;
		double oldLL = theta.LL;
		do{
			oldTheta = theta.theta;
			_NROnlyTheta();
			++i;
		} while(i<numThetaOnlyUpdates && theta.theta != oldTheta);

		//update all params
		_NRAllParams();

		//e) do we break EM? Check LL
		theta.LL = data->calcLogLikelihood(pGenotype);
		if((theta.LL - oldLL) < maxEpsilon)
			break;

		//maybe theta = 0?
		if(theta.theta < 0.1/(double) data->size()){
			//(theta is somewhere between 1/numLoci and 0,
			oldTheta = theta.theta;

			//test with theta = 0.0
			theta.setTheta(0.0);
			fillPGenotype(pGenotype, theta);
			theta.LL = data->calcLogLikelihood(pGenotype);

			//theta is between zero and a very small number -> don't care about exact value
			if(theta.LL < oldLL){
				theta.setTheta(oldTheta);
				theta.LL = oldLL;
			}
			break;
		}

		if(extraVerbose)
			logfile->list(toString(iter) + ") current theta = " + toString(theta.theta));
		//For debugging
		//std::cout << std::setprecision(9) << iter << ") theta = " << thetaContainer.theta << "\tLL = " << thetaContainer.LL << "\teps = " << fabs(oldLL - thetaContainer.LL) << std::endl;
	}
	if(extraVerbose) logfile->list("EM converged, current theta = " + toString(theta.theta));

}

void TThetaEstimator::_estimateConfidenceInterval(){
	using namespace genometools;
	//we estimate an approximate confidence interval for theta using the Fisher information
	//This function assumes that EM has already been run!

	//calculate P(g|theta, pi)
	fillPGenotype(pGenotype, theta);

	//calclate d/dtheta P(g|theta, pi)
	TGenotypeData deriv_pGenotype;

	for(Base k = Base::min; k < Base::max; ++k){
		//homozygous genotype
		const auto hom = genotype(k, k);
		deriv_pGenotype[hom] = (theta.baseFreq[k].get() * theta.baseFreq[k].get() - theta.baseFreq[k].get()) * theta.expTheta;
		//heterozygous genotypes
		for(Base l = next(k); l < Base::max; ++l){
			const auto het = genotype(k, l);
			deriv_pGenotype[het] = 2.0 * theta.baseFreq[k].get() * theta.baseFreq[l].get() * theta.expTheta;
		}
	}

	double FisherInfo = _calcFisherInfo(pGenotype, deriv_pGenotype);

	//estimate confidence interval
	//TODO: Fisher Info can be negative -> SQRT will be nan!
	theta.thetaConfidence = 1.96 / sqrt(FisherInfo);
}

//------------------------------------------------------------
//Functions to run estimation-
//------------------------------------------------------------
bool TThetaEstimator::estimateTheta(){
	estimationSuccessful = false;
	if(data->sizeWithData() < minSitesWithData){
		logfile->warning("Can not estimate theta: only " + toString(data->sizeWithData()) + " sites with data in this region (minSitesWithData = " + toString(minSitesWithData) + ")!");
		return false;
	}

	//estimate starting parameters
	findGoodStartingTheta(data, theta, "");

	//Run EM
	if(extraVerbose){
		logfile->startIndent("Running EM to find ML estimate:");
		_runEMForTheta();
		logfile->endIndent();
	} else {
		logfile->listFlush("Running EM to find ML estimate ...");
		_runEMForTheta();
		logfile->done();
	}
	logfile->conclude("theta was estimated at ", theta.theta);

	//confidence intervals
	logfile->listFlush("Estimating approximate confidence intervals from Fisher-Information ...");
	_estimateConfidenceInterval();
	logfile->done();
	logfile->conclude("95% confidence intervals are theta +- " + toString(theta.thetaConfidence));
	estimationSuccessful = true;
	return true;
}

void TThetaEstimator::setTheta(const double Theta){
	theta.setTheta(Theta);
};

void TThetaEstimator::setBaseFreq(const GenotypeLikelihoods::TBaseProbabilities & BaseFreq){
	theta.baseFreq = BaseFreq;
};

void TThetaEstimator::addToHeader(std::vector<std::string> & header, std::string prefix){
	data->addToHeader(header, prefix);
	header.push_back(prefix + "pi(A)");
	header.push_back(prefix + "pi(C)");
	header.push_back(prefix + "pi(G)");
	header.push_back(prefix + "pi(T)");
	header.push_back(prefix + "theta_MLE");
	header.push_back(prefix + "theta_C95_l");
	header.push_back(prefix + "theta_C95_u");
	header.push_back(prefix + "LL");
}

void TThetaEstimator::writeEstimateFrequenciesAndTheta(coretools::TOutputFile & out){
	if(estimationSuccessful){
		//base frequencies
		for(genometools::Base k = genometools::Base::min; k < genometools::Base::max; ++k){
			out << theta.baseFreq[k];
		}

		//theta estimates
		out << theta.theta;
		out	<< theta.theta - theta.thetaConfidence;
		out	<< theta.theta + theta.thetaConfidence;
		out	<< theta.LL;
	} else {
		//frequencies
		out << "-" << "-" << "-" << "-";

		//theta estimates
		out << "-" << "-" << "-" << "-";
	}
};

void TThetaEstimator::writeResultsToFile(coretools::TOutputFile & out){
	//number of sites
	data->writeSite(out);

	//estimated params
	writeEstimateFrequenciesAndTheta(out);
};

void TThetaEstimator::calcLikelihoodSurface(coretools::TOutputFile & out, uint32_t steps){
	//write header
	out.writeHeader({"log10(theta)", "theta", "LL"});
	out.setPrecision(12);

	//calculate likelihood surface
	double minLogTheta = -5.0;
	double maxLogTheta = -1.0;
	double stepSize = (maxLogTheta - minLogTheta) / ((double) steps - 1.0);

	for(uint32_t i=0; i<steps; ++i){
		//calc theta and expTheta
		theta.setLogTheta(minLogTheta + stepSize*i);

		//calculate	substitution probabilities and Likelihood
		fillPGenotype(pGenotype, theta);
		theta.LL = data->calcLogLikelihood(pGenotype);

		//write results
		out << theta.logTheta << theta.theta << theta.LL << std::endl;
	}
};

void TThetaEstimator::bootstrapTheta(){
	logfile->listFlush("Bootstrapping sites ...");

	data->bootstrap(*randomGenerator);
	logfile->done();

	//estimate theta
	estimateTheta();

	//clean up
	data->clearBootstrap();
};

//---------------------------------------------------------------
//TThetaEstimatorRatio
//---------------------------------------------------------------
TThetaEstimatorRatio::TThetaEstimatorRatio(coretools::TParameters & params, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator):TThetaEstimator_base(params, Logfile, RandomGenerator){
	clearCounters();

	//data2
	if(useTmpFile){
		data2 = new TThetaEstimatorDataFile(tmpFileName + "2.tmp.gz");
	} else
		data2 = new TThetaEstimatorDataVector();
	data2Initialized = true;

	//MCMC params
	logfile->startIndent("Parameters of MCMC algorithm:");
	burnin = params.getParameterWithDefault<int>("burnin", 1000);
	logfile->list("Will run a burnin of " + toString(burnin) + " iterations.");
	numIterations = params.getParameterWithDefault<int>("iterations", 10000);
	logfile->list("Will run MCMC for " + toString(numIterations) + " iterations.");
	thinning = params.getParameterWithDefault<int>("thinning", 1);
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

	//normal prior on ratio phi = log(theta_1 / theta_2)
	phiPriorMean = params.getParameterWithDefault("phiPriorMean", 0.0);
	phiPriorVar = params.getParameterWithDefault("phiPriorVar", 1.0);
	phiPriorOneOverTwoVar = 1.0 / 2.0 / phiPriorVar;
	logfile->list("Will assume a normal prior on phi ~ N(" + toString(phiPriorMean) + ", " + toString(phiPriorVar) + ").");

	//proposal kernel
	sdProposalKernelTheta1 = params.getParameterWithDefault("sdProposalTheta", 0.1);
	sdProposalKernelTheta2 = sdProposalKernelTheta1;
	sdProposalKernelBaseFreq1 = params.getParameterWithDefault("sdProposalFreq", 0.01);
	sdProposalKernelBaseFreq2 = sdProposalKernelBaseFreq1;
	logfile->list("Will use initial proposal kernel standard deviations of " + toString(sdProposalKernelTheta1) + " and " + toString(sdProposalKernelBaseFreq1) + " for thetas and base frequencies, respectively.");
	logfile->endIndent();

	//params regarding initial search
	readParametersRegardingInitialSearch(params);
};

void TThetaEstimatorRatio::clearCounters(){
	numAcceptedTheta1 = 0;
	numAcceptedTheta2 = 0;
	numAcceptedBaseFreq1 = 0;
	numAcceptedBaseFreq2 = 0;
};

void TThetaEstimatorRatio::concludeAcceptanceRate(int numAccepted, int length, std::string name){
	double acceptanceRate = (double) numAccepted / (double) length;
	logfile->conclude("Acceptance rate " + name + " = " + toString(acceptanceRate));
};

void TThetaEstimatorRatio::concludeAcceptanceRateUpdateProposal(int numAccepted, int length, double & sd, std::string name){
	double acceptanceRate = (double) numAccepted / (double) length;
	sd *= acceptanceRate * 3.0;
	logfile->conclude("Acceptance rate " + name + " = " + toString(acceptanceRate) + " (updated proposal sd to " + toString(sd) + ")");
};

void TThetaEstimatorRatio::concludeAcceptanceRates(int length){
	concludeAcceptanceRate(numAcceptedTheta1, length, "theta 1");
	concludeAcceptanceRate(numAcceptedTheta2, length, "theta 2");
	concludeAcceptanceRate(numAcceptedBaseFreq1, length, "base frequencies 1");
	concludeAcceptanceRate(numAcceptedBaseFreq2, length, "base frequencies 1");
};

void TThetaEstimatorRatio::concludeAcceptanceRatesUpdateProposal(int length){
	concludeAcceptanceRateUpdateProposal(numAcceptedTheta1, length, sdProposalKernelTheta1, "theta 1");
	concludeAcceptanceRateUpdateProposal(numAcceptedTheta2, length, sdProposalKernelTheta2, "theta 2");
	concludeAcceptanceRateUpdateProposal(numAcceptedBaseFreq1, length, sdProposalKernelBaseFreq1, "base frequencies 1");
	concludeAcceptanceRateUpdateProposal(numAcceptedBaseFreq2, length, sdProposalKernelBaseFreq2, "base frequencies 1");
};

void TThetaEstimatorRatio::estimateRatio(std::string ouputName){
	logfile->startIndent("Running MCMC to estimate phi = log(theta1 / theta2):");

	//check if there is sufficient data
	logfile->list(toString(data->sizeWithData()) + " sites with data available for region 1.");
	if(data->numSitesWithData < minSitesWithData)
		throw "Not enough sites for region 1!";
	logfile->list(toString(data2->sizeWithData()) + " sites with data available for region 2.");
	if(data2->numSitesWithData < minSitesWithData)
		throw "Not enough sites for region 2!";

	//get good starting values
	findGoodStartingTheta(data, theta, " region 1");
	findGoodStartingTheta(data2, theta2, " region 2");

	//first run burnin
	clearCounters();
	if(burnin > 0){
		int oldProg = 0;
		std::string progressString = "Running burnin of length " + toString(burnin) + " ...";
		logfile->listFlush(progressString + "(0%)");
		for(int i=0; i<burnin; ++i){
			oneMCMCIteration();

			//print progress
			int prog = (double) i / (double) burnin * 100.0;
			if(prog > oldProg){
				oldProg = prog;
				logfile->overListFlush(progressString, "(", oldProg, "%)");
			}
		}
		logfile->overList(progressString + "done!   ");
		concludeAcceptanceRatesUpdateProposal(burnin);
	}


	//open MCMC output file
	ouputName += "_thetaRatioMCMC.txt.gz";
	logfile->list("Will write MCMC chain to file '" + ouputName + "'.");
	gz::ogzstream out(ouputName.c_str());
	if(!out) throw "Failed to open file '" + ouputName + "' for writing!";

	//write header
	out << "log_theta_1\tlog_theta_2\tlog_phi\n";

	//now run chain with sampling
	clearCounters();
	int oldProg = 0;
	std::string progressString = "Running MCMC chain of length " + toString(numIterations) + " ...";
	logfile->listFlush(progressString + "(0%)");
	for(int i=0; i<numIterations; ++i){
		oneMCMCIteration();

		//print to file
		if(i % thinning == 0){
			out << theta.logTheta << "\t" << theta2.logTheta << "\t" << theta.logTheta - theta2.logTheta << "\n";
		}

		//print progress
		int prog = (double) i / (double) numIterations * 100.0;
		if(prog > oldProg){
			oldProg = prog;
			logfile->overListFlush(progressString, "(", oldProg, "%)");
		}
	}
	logfile->overList(progressString + " done!   ");
	concludeAcceptanceRates(numIterations);

	//clean up
	logfile->endIndent();
	out.close();
};

bool TThetaEstimatorRatio::updateTheta(TThetaEstimatorData* thisData, Theta & thisTheta, double otherLogThetaMean, double thisSdProposalKernel){
	//propose
	double newLogTheta = randomGenerator->getNormalRandom(thisTheta.logTheta, thisSdProposalKernel);
	double newExpTheta = exp(-exp(newLogTheta)); //we update log(theta) but need exp(-theta)

	//calc LL
	fillPGenotype(pGenotype, newExpTheta, thisTheta.baseFreq);
	double newLL = thisData->calcLogLikelihood(pGenotype);

	//calc hastings ratio with prior
	//we use a uniform prior on log(theta)
	//and normal on log(phi) = log(theta) - log(theta2)
	double logH = newLL - thisTheta.LL ;
	double newLogPhiMinusMean = newLogTheta - otherLogThetaMean;
	double oldLogPhiMinusMean = thisTheta.logTheta - otherLogThetaMean;
	logH += phiPriorOneOverTwoVar * (newLogPhiMinusMean*newLogPhiMinusMean - oldLogPhiMinusMean*oldLogPhiMinusMean);

	//accept or reject
	if(log(randomGenerator->getRand()) < logH){
		thisTheta.setLogTheta(newLogTheta, newLL);
		return true;
	} else return false;
};

bool TThetaEstimatorRatio::updateBaseFrequencies(TThetaEstimatorData* thisData, Theta & thisTheta, double thisSdProposalKernel){
	using namespace genometools;
	//propose: select one frequency at random and shift this one
	//make sure frequencies are not outside [0,1]
	const auto thisBase = Base(randomGenerator->sample(4));
	TBaseProbabilities tmpBaseFreq;
	double tmp = thisTheta.baseFreq[thisBase].get() + randomGenerator->getNormalRandom(0.0, thisSdProposalKernel);
	if(tmp > 1.0){
		tmpBaseFreq[thisBase] = 2.0 - tmp;
	} else if (tmp < 0.0){
		tmpBaseFreq[thisBase] = -tmp;
	} else {
		tmpBaseFreq[thisBase] = tmp;
	}

	//now scale all others so the sum will be 1.0
	double scale = (double) tmpBaseFreq[thisBase].complement() / (double) thisTheta.baseFreq[thisBase].complement();
	for(Base k = Base::min; k < Base::max; ++k){
		if(k != thisBase){
			tmpBaseFreq[k] = thisTheta.baseFreq[thisBase].get() * scale;
		}
	}

	//calc LL & hastings ratio (use uniform prior, i.e. all combinations are equally likely)
	fillPGenotype(pGenotype, thisTheta.expTheta, tmpBaseFreq);
	double newLL = thisData->calcLogLikelihood(pGenotype);
	double logH = newLL - thisTheta.LL;

	//accept or reject
	if(log(randomGenerator->getRand()) < logH){
		thisTheta.baseFreq = tmpBaseFreq;
		thisTheta.LL = newLL;
		return true;
	} else return false;
};

void TThetaEstimatorRatio::oneMCMCIteration(){
	//update theta
	numAcceptedTheta1 += updateTheta(data, theta, theta2.logTheta + phiPriorMean, sdProposalKernelTheta1);
	numAcceptedTheta2 += updateTheta(data2, theta2, theta.logTheta - phiPriorMean, sdProposalKernelTheta2);

	//update base frequencies
	numAcceptedBaseFreq1 += updateBaseFrequencies(data, theta, sdProposalKernelBaseFreq1);
	numAcceptedBaseFreq2 += updateBaseFrequencies(data2, theta2, sdProposalKernelBaseFreq2);
};

}; //end namespace
