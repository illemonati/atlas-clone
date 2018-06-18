/*
 * TThetaEstimator.cpp
 *
 *  Created on: Sep 25, 2017
 *      Author: phaentu
 */

#include "TThetaEstimator.h"



//---------------------------------------------------------------
//TThetaEstimator_base
//---------------------------------------------------------------
TThetaEstimator_base::TThetaEstimator_base(TLog* Logfile){
	logfile = Logfile;
	numGenotypes = 10;
	initTmpStorage();

	useTmpFile = false;
	minSitesWithData = 0;
	data = NULL;
	dataInitialized = false;
	extraVerbose = false;
};

TThetaEstimator_base::TThetaEstimator_base(TParameters & params, TLog* Logfile){
	logfile = Logfile;
	numGenotypes = 10;
	initTmpStorage();

	//data
	useTmpFile = params.parameterExists("useTmpFile");
	if(useTmpFile){
		tmpFileName = params.getParameterStringWithDefault("useTmpFile", "temporaryDataForThetaEstimation");
		logfile->list("Will write temporary data to file(s) with prefix '" + tmpFileName + "'.");
		data = new TThetaEstimatorDataFile(numGenotypes, tmpFileName + ".tmp.gz");
	} else
		data = new TThetaEstimatorDataVector(numGenotypes);
	dataInitialized = true;

	//minimum window size
	minSitesWithData = params.getParameterIntWithDefault("minSitesWithData", 1000);

	//extra verbosity
	extraVerbose = params.parameterExists("extraVerbose");
};

void TThetaEstimator_base::initTmpStorage(){
	//initialize arrays
	pGenotype = new double[numGenotypes];
	baseFreq = new double[4];
}

void TThetaEstimator_base::readParametersRegardingInitialSearch(TParameters & params){
	initalTheta = params.getParameterDoubleWithDefault("initTheta", 0.01);
	logfile->list("Will start with an initial theta of " + toString(initalTheta) + ".");
	initThetaNumSearchIterations = params.getParameterDoubleWithDefault("initThetaNumSearchIterations", 10);
	if(initThetaNumSearchIterations > 0){
		logfile->list("Will run " + toString(initThetaNumSearchIterations) + " iterations of a crude search for an initial theta.");
		initThetaSearchFactor = params.getParameterDoubleWithDefault("initThetaSearchFactor", 100);
		logfile->list("The initial search factor will be " + toString(initThetaSearchFactor) + ".");
	} else {
		initThetaSearchFactor = 0;
	}
};


void TThetaEstimator_base::fillPGenotype(double* & pGeno, double & expTheta, double* baseFrequencies){
	//assumes that base frequencies are set!
	static int i;
	static int j;
	for(i=0; i<4; ++i){
		//homozygous genotypes
		pGeno[genoMap.getGenotype(i,i)] = baseFrequencies[i] * (expTheta + baseFrequencies[i] * (1.0 - expTheta));
		//heterozygous genotypes
		for(j=i+1; j<4; ++j){
			pGeno[genoMap.getGenotype(i,j)] = 2.0 * baseFrequencies[i] * baseFrequencies[j] *  (1.0 - expTheta);
		}
	}
};

void TThetaEstimator_base::findGoodStartingTheta(TThetaEstimatorData* thisData, Theta & thisTheta, std::string tag){
	logfile->listFlush("Estimating initial parameters" + tag + " ...");

	//set base frequencies to initial base frequencies
	thisData->fillBaseFreq(baseFreq);

	//variables
	double initTheta = initalTheta;
	double oldTheta = initTheta;
	double expTheta = exp(-initTheta);

	//calc initial LL
	fillPGenotype(pGenotype, expTheta, baseFreq);
	thisTheta.LL = data->calcLogLikelihood(pGenotype);

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
			fillPGenotype(pGenotype, expTheta, baseFreq);
			thisTheta.LL = data->calcLogLikelihood(pGenotype);
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
				fillPGenotype(pGenotype, expTheta, baseFreq);
				thisTheta.LL = data->calcLogLikelihood(pGenotype);
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
	logfile->conclude("Initial base frequencies: Pi(A) = " + toString(baseFreq[0]) + ", Pi(C) = " + toString(baseFreq[1]) + ", Pi(G) = " + toString(baseFreq[2]) + ", Pi(T) = " + toString(baseFreq[3]));
	logfile->conclude("Initial theta = ", theta.theta);
}

//-------------------------------------------------------
//TThetaEstimator
//-------------------------------------------------------
TThetaEstimator::TThetaEstimator(TParameters & params, TLog* Logfile):TThetaEstimator_base(params, Logfile){
	initAdditionalTmpStorage();

	//parse
	logfile->startIndent("Parameters of EM algorithm:");
	numIterations = params.getParameterIntWithDefault("iterations", 100);
	logfile->list("Will run up to " + toString(numIterations) + " iterations.");
	numThetaOnlyUpdates = params.getParameterIntWithDefault("iterationsThetaOnly", 10);
	logfile->list("In each iteration, theta will be updated " + toString(numThetaOnlyUpdates) + " times.");

	maxEpsilon = params.getParameterDoubleWithDefault("maxEps", 0.000001);
	logfile->list("Will run EM until deltaLL < " + toString(maxEpsilon) + ".");
	NewtonRaphsonNumIterations = params.getParameterIntWithDefault("NRiterations", 10);
	logfile->list("Will run Newton-Raphson algorithm up to " + toString(NewtonRaphsonNumIterations) + " times.");
	NewtonRaphsonMaxF = params.getParameterDoubleWithDefault("maxF", 0.00001);
	logfile->list("Will run Newton-Raphson algorithm until max(F) < " + toString(NewtonRaphsonMaxF) + ".");

	//params regarding initial search
	readParametersRegardingInitialSearch(params);

	//done
	logfile->endIndent();
}

TThetaEstimator::TThetaEstimator(TLog* Logfile):TThetaEstimator_base(logfile){
	initAdditionalTmpStorage();

	//set EM parameters to default
	numIterations = -1;
	numThetaOnlyUpdates = -1;
	maxEpsilon = 0.0;
	NewtonRaphsonNumIterations = -1;
	NewtonRaphsonMaxF = 0.0;
	initalTheta = 0.0;
	initThetaSearchFactor = -1;
	initThetaNumSearchIterations = -1;
};

void TThetaEstimator::initAdditionalTmpStorage(){
	//initialize arrays
	P_G = new double[numGenotypes];
	P_g_oneSite = new double[numGenotypes];
};

void TThetaEstimator::clear(){
	data->clear();
};

void TThetaEstimator::add(TSite & site){
	data->add(site);
};

void TThetaEstimator::fillP_G(){
	//assumes that pGenotype is set!
	for(int g=0; g<numGenotypes; ++g)
		P_G[g] = 0.0;

	//calculate P_g for each site
	double* d;
	data->begin();
	do{
		double sum = 0.0;
		d = data->curGenotypeLikelihoods();
		for(int g=0; g<numGenotypes; ++g){
			P_g_oneSite[g] =  d[g] * pGenotype[g];
			sum += P_g_oneSite[g];
		}
		for(int g=0; g<numGenotypes; ++g)
			P_G[g] += P_g_oneSite[g] / sum;

	} while(data->next());
};

double TThetaEstimator::calcFisherInfo(double* _pGenotype, double* deriv_pGenotype){
//sum Ri over all sites
	double FisherInfo = 0.0;
	double Ri, Ri_a, Ri_b;
	double* d;

	data->begin();
	do{
		//calc Ri
		Ri_a = 0.0; Ri_b = 0.0;
		d = data->curGenotypeLikelihoods();
		for(int g=0; g<numGenotypes; ++g){
			Ri_a += d[g] * deriv_pGenotype[g];
			Ri_b += d[g] * _pGenotype[g];
		}
		Ri = Ri_a / Ri_b;

		//add to Fisher Info
		FisherInfo += Ri * (Ri + 1.0);
	} while(data->next());

	return(FisherInfo);
};


void TThetaEstimator::runEMForTheta(){
	//prepare storage
	double tmp[4];
	double tmpSum;
	arma::mat Jacobian(6,6);
	arma::vec F(6);
	arma::mat JxF(6,6);
	Genotype geno;
	double maxF;
	int failedAttempts = 0;
	double oldTheta, rho, mu = data->sizeWithData();
	double oldLL = -9e100;

	//start EM loop
	int numThetaOnlyUpdatesDone = numThetaOnlyUpdates; //do regular step first
	numThetaOnlyUpdatesDone = 0;
	int totIterations = numIterations * numThetaOnlyUpdates;
	for(int iter = 0; iter < totIterations; ++iter){
		//a) pre-calc expTheta
		oldTheta = theta.theta;
		rho = theta.expTheta / (1.0 - theta.expTheta);

		//b) calculate	substitution probabilities
		fillPGenotype(pGenotype, theta.expTheta, baseFreq);

		//c) Calculate all genotype probabilities for all sites
		fillP_G();
		//data->fillP_G(P_G, pGenotype);

		//d) Find new parameter estimates using Newton-Raphson
		if(numThetaOnlyUpdatesDone < numThetaOnlyUpdates){
			//update only theta: most difficult parameter and it is much faster to update only this one alone.
			for(int n=0; n<NewtonRaphsonNumIterations; ++n){
				//i) calculate F() (Note: index is zero based!)
				F(4) = data->sizeWithData();
				for(int k=0; k<4; ++k){
					geno = genoMap.getGenotype(k, k);
					F(4) -= P_G[geno] * (rho + 1.0 ) / (rho + baseFreq[k]);
				}
				//ii) fill Jacobian (Note: index is zero based!)
				Jacobian(4,4) = 0.0;
				for(int k=0; k<4; ++k){
					tmpSum = P_G[genoMap.getGenotype(k, k)] / ((baseFreq[k] + rho)*(baseFreq[k] + rho));
					Jacobian(4,4) += tmpSum * (1.0 - baseFreq[k]);
				}

				//iii) now estimate new parameters
				rho = rho - F(4) / Jacobian(4,4);

				//check if we break
				if(F(4) < NewtonRaphsonMaxF){
					theta.setTheta(-log(rho / (1.0 + rho)));
					break;
				}
			}
			++numThetaOnlyUpdatesDone;
			if(theta.theta == oldTheta) numThetaOnlyUpdatesDone = numThetaOnlyUpdates;
		} else {
			numThetaOnlyUpdatesDone = 0;
			//update all parameters in EM
			for(int n=0; n<NewtonRaphsonNumIterations; ++n){
				//i) calculate F (Note: index is zero based!)
				F(4) = data->sizeWithData();
				F(5) = 0.0;
				for(int k=0; k<4; ++k){
					geno = genoMap.getGenotype(k, k);
					tmpSum = 0.0;
					for(int l=0; l<4; ++l){
						if(l != k){
							tmpSum += P_G[genoMap.getGenotype(k, l)];
						}
					}
					F(k) = P_G[geno] * (1.0 + baseFreq[k] / (rho + baseFreq[k])) + tmpSum - mu * baseFreq[k];
					F(4) -= P_G[geno] * (rho + 1.0 ) / (rho + baseFreq[k]);
					F(5) += baseFreq[k];
				}
				F(5) = F(5) - 1.0;

				//ii) fill Jacobian (Note: index is zero based!)
				Jacobian.zeros();
				tmpSum = 0.0;
				for(int k=0; k<4; ++k){
					tmp[k] = P_G[genoMap.getGenotype(k, k)] / ((baseFreq[k] + rho)*(baseFreq[k] + rho));
					tmpSum += tmp[k];
				}

				for(int k=0; k<4; ++k){
					Jacobian(k,k) = tmp[k] * rho - mu;
					Jacobian(k,4) = - tmp[k];
					Jacobian(5,k) = 1.0;
					Jacobian(4,k) = tmp[k] * (rho + 1.0);
					Jacobian(k,5) = - baseFreq[k];
					Jacobian(4,4) += tmp[k] * (1.0 - baseFreq[k]);
				}

				//iii) now estimate new parameters
				if(solve(JxF, Jacobian, F)){
					for(int k=0; k<4; ++k){
						baseFreq[k] -= JxF(k);
					}
					rho -= JxF(4);
					mu -= JxF(5);

					//check if we break
					maxF = 0.0;
					for(int i=0; i<6; ++i){
						if(F(i) > maxF) maxF = F(i);
					}

					if(maxF < NewtonRaphsonMaxF || n == (NewtonRaphsonNumIterations-1)){
						theta.setTheta(-log(rho / (1.0 + rho)));
						break;
					}
				} else {
					++failedAttempts;

					//solve did not work -> start with higher theta!
					theta.setTheta(initalTheta);
					for(int i=0; i<failedAttempts; ++i)
						theta.theta *= 10.0;

					//reset others
					mu = data->sizeWithData();
					theta.LL = -9e100;
					iter = 0;
					numThetaOnlyUpdatesDone = 0;
					break;
				}
			}
		}

		//e) do we break EM? Check LL
		if(iter > 0 && iter % numThetaOnlyUpdates == 0){
			oldLL = theta.LL;
			theta.LL = data->calcLogLikelihood(pGenotype);
			if(theta.LL > -9e100 && (theta.LL - oldLL) < maxEpsilon) break;

			//maybe theta = 0?
			if(theta.theta < 0.1/(double) data->size()){
				oldLL = theta.LL;
				oldTheta = theta.theta;

				//test with theta = 0.0
				theta.setTheta(0.0);
				fillPGenotype(pGenotype, theta.expTheta, baseFreq);
				theta.LL = data->calcLogLikelihood(pGenotype);

				if(theta.LL < oldLL){
					theta.setTheta(oldTheta);
					theta.LL = oldLL;
				}
				break;
			}
		}

		if(extraVerbose) logfile->list(toString(iter) + ") current theta = " + toString(theta.theta));
		//For debugging
		//std::cout << std::setprecision(9) << iter << ") theta = " << thetaContainer.theta << "\tLL = " << thetaContainer.LL << "\teps = " << fabs(oldLL - thetaContainer.LL) << std::endl;
	}
	if(extraVerbose) logfile->list("EM converged, current theta = " + toString(theta.theta));

}

void TThetaEstimator::estimateConfidenceInterval(){
	//we estimate an approximate confidence interval for theta using the Fisher information
	//This function assumes that EM has already been run!

	//calculate P(g|theta, pi)
	fillPGenotype(pGenotype, theta.expTheta, baseFreq);

	//calclate d/dtheta P(g|theta, pi)
	double deriv_pGenotype[numGenotypes];
	for(int k=0; k<4; ++k){
		//homozygous genotype
		deriv_pGenotype[genoMap.getGenotype(k, k)] = (baseFreq[k] * baseFreq[k] - baseFreq[k]) * theta.expTheta;
		//heterozygous genotypes
		for(int l=k+1; l<4; ++l){
			deriv_pGenotype[genoMap.getGenotype(k, l)] = 2.0 * baseFreq[k] * baseFreq[l] * theta.expTheta;
		}
	}

	double FisherInfo = calcFisherInfo(pGenotype, deriv_pGenotype);

	//estimate confidence interval
	//TODO: Fisher Info can be negative -> SQRT will be nan!
	theta.thetaConfidence = 1.96 / sqrt(FisherInfo);
}

//------------------------------------------------------------
//Functions to run estimation-
//------------------------------------------------------------
bool TThetaEstimator::estimateTheta(){
	if(data->sizeWithData() < minSitesWithData){
		logfile->write("Can not estimate theta, less than minSitesWithData = " + toString(minSitesWithData) + " sites with data in this region!");
		return false;
	}

	//estimate starting parameters
	findGoodStartingTheta(data, theta, "");

	//Run EM
	if(extraVerbose){
		logfile->startIndent("Running EM to find ML estimate:");
		runEMForTheta();
		logfile->endIndent();
	} else {
		logfile->listFlush("Running EM to find ML estimate ...");
		runEMForTheta();
		logfile->done();
	}
	logfile->conclude("theta was estimated at ", theta.theta);

	//confidence intervals
	logfile->listFlush("Estimating approximate confidence intervals from Fisher-Information ...");
	estimateConfidenceInterval();
	logfile->done();
	logfile->conclude("95% confidence intervals are theta +- " + toString(theta.thetaConfidence));
	return true;
}

void TThetaEstimator::setTheta(double Theta){
	theta.setTheta(Theta);
}

void TThetaEstimator::setBaseFreq(TBaseFrequencies & BaseFreq){
	for(int i=0; i<4; ++i)
		baseFreq[i] = BaseFreq[i];
}

void TThetaEstimator::writeHeader(std::ofstream & out){
	data->writeHeader(out);
	out << "\tpi(A)\tpi(C)\tpi(G)\tpi(T)\ttheta_MLE\ttheta_C95_l\ttheta_C95_u\tLL";
}

void TThetaEstimator::writeThetas(std::ofstream & out){
	out << "\t" << theta.theta;
	out	<< "\t" << theta.theta - theta.thetaConfidence;
	out	<< "\t" << theta.theta + theta.thetaConfidence;
	out	<< "\t" << theta.LL;
	out	<< std::endl;
}
void TThetaEstimator::writeResultsToFile(std::ofstream & out){
	//number of sites
	data->writeSize(out);

	//estimated params
	for(int i=0; i<4; ++i)
		out << "\t" << baseFreq[i];
	writeThetas(out);
}

void TThetaEstimator::calcLikelihoodSurface(std::ofstream & out, int & steps){
	//write header
	out << "log10(theta)\ttheta\tLL\n";

	//calculate likelihood surface
	double minLogTheta = -5.0;
	double maxLogTheta = -1.0;
	double stepSize = (maxLogTheta - minLogTheta) / ((double) steps - 1.0);
	double theta;
	double LL;
	double expTheta;
	double logTheta;

	for(int i=0; i<steps; ++i){
		//calc theta and expTheta
		logTheta = minLogTheta + stepSize*i;
		theta = pow(10.0, logTheta);
		expTheta = exp(-theta);

		//calculate	substitution probabilities and Likelihood
		fillPGenotype(pGenotype, expTheta, baseFreq);
		LL = data->calcLogLikelihood(pGenotype);

		//write results
		out << std::setprecision(12) << logTheta << "\t" << theta << "\t" << LL << "\n";
	}
}

void TThetaEstimator::bootstrapTheta(TRandomGenerator & randomGenerator, std::ofstream & out){
	logfile->listFlush("Bootstrapping sites ...");

	data->bootstrap(randomGenerator);
	logfile->done();

	//estimate theta
	estimateTheta();

	//write results
	writeResultsToFile(out);

	//clean up
	data->clearBootstrap();
}


//---------------------------------------------------------------
//TThetaEstimatorRatio
//---------------------------------------------------------------
TThetaEstimatorRatio::TThetaEstimatorRatio(TParameters & params, TLog* Logfile):TThetaEstimator_base(params, Logfile){
	initAdditionalTmpStorage();

	//data2
	if(useTmpFile){
		data2 = new TThetaEstimatorDataFile(numGenotypes, tmpFileName + "2.tmp.gz");
	} else
		data2 = new TThetaEstimatorDataVector(numGenotypes);
	data2Initialized = true;

	//MCMC params
	logfile->startIndent("Parameters of MCMC algorithm:");
	numIterations = params.getParameterIntWithDefault("iterations", 10000);
	logfile->list("Will run MCMC for " + toString(numIterations) + " iterations.");
	burnin = params.getParameterIntWithDefault("burnin", 1000);
	logfile->list("Will run a burnin of " + toString(burnin) + " iterations.");

	//normal prior on ratio phi = log(theta_1 / theta_2)
	phiPriorMean = params.getParameterDoubleWithDefault("phiPriorMean", 0.0);
	phiPriorVar = params.getParameterDoubleWithDefault("phiPriorVar", 1.0);
	phiPriorOneOverTwoVar = 1.0 / 2.0 / phiPriorVar;
	logfile->list("Will assume a normal prior on phi ~ N(" + toString(phiPriorMean) + ", " + toString(phiPriorVar) + ").");

	//proposal kernel
	sdProposalKernel = params.getParameterDoubleWithDefault("sdProposal", 0.1);

	//params regarding initial search
	readParametersRegardingInitialSearch(params);

	//done
	logfile->endIndent();
};

void TThetaEstimatorRatio::initAdditionalTmpStorage(){
	baseFreq2 = new double[4];
	tmpBaseFreq = new double[4];
};

void TThetaEstimatorRatio::estimateRatio(TRandomGenerator & randomGenerator, std::string ouputName){
	logfile->startIndent("Running MCMC to estimate phi = log(theta1 / theta2):");

	//open MCMC output file
	ouputName += "_thetaRatioMCMC.txt.gz";
	logfile->list("Will write MCMC chain to file '" + ouputName + "'.");
	gz::ogzstream out(ouputName.c_str());
	if(!out) throw "Failed to open file '" + ouputName + "' for writing!";

	//get good starting values
	findGoodStartingTheta(data, theta, " region 1");
	findGoodStartingTheta(data2, theta2, " region 2");

	//now run MCMC





	//clean up
	logfile->endIndent();
	out.close();


};

void TThetaEstimatorRatio::updateTheta(TThetaEstimatorData* thisData, Theta & thisTheta, double* theseBaseFrequencies, double otherLogThetaMean, TRandomGenerator & randomGenerator){
	//propose
	double newLogTheta = randomGenerator.getNormalRandom(thisTheta.logTheta, sdProposalKernel);
	double newExpTheta = exp(-exp(newLogTheta)); //we update log(theta) but need exp(-theta)

	//calc LL
	fillPGenotype(pGenotype, newLogTheta, theseBaseFrequencies);
	double newLL = thisData->calcLogLikelihood(pGenotype);

	//calc hastings ratio with prior
	//we use a uniform prior on log(theta)
	//and normal on log(phi) = log(theta) - log(theta2)
	double logH = newLL - thisTheta.LL ;
	double newLogPhiMinusMean = newLogTheta - otherLogThetaMean;
	double oldLogPhiMinusMean = thisTheta.logTheta - otherLogThetaMean;
	logH += phiPriorOneOverTwoVar * (oldLogPhiMinusMean*oldLogPhiMinusMean - newLogPhiMinusMean*newLogPhiMinusMean);

	//accept or reject
	if(log(randomGenerator.getRand()) < logH)
		thisTheta.setLogTheta(newLogTheta, newLL);
};

void TThetaEstimatorRatio::updateBaseFrequencies(TThetaEstimatorData* thisData, Theta & thisTheta, double* theseFrequencies, TRandomGenerator & randomGenerator){
	//propose: select one frequency at random and shift this one
	//make sure frequencies are not outsiue [0,1]
	int numOutsideRange = 1;
	while(numOutsideRange > 0){
		int thisFreq = randomGenerator.pickOne(4);
		double delta = randomGenerator.getRand() * 0.01 - 0.005; //shift at max by +/- 0.005

		tmpBaseFreq[thisFreq] = theseFrequencies[thisFreq] + delta;

		numOutsideRange = 0;

		for(int i=0; i<4; ++i){
			if(i != thisFreq)
				tmpBaseFreq[i] = theseFrequencies[i] - delta / 3.0;
			if(tmpBaseFreq[i] < 0.0 || tmpBaseFreq[i] > 1.0)
				++numOutsideRange;
		}
	}

	//calc LL & hastings ratio (use uniform prior, i.e. all combinations are equally likely)
	fillPGenotype(pGenotype, thisTheta.expTheta, tmpBaseFreq);
	double newLL = thisData->calcLogLikelihood(pGenotype);
	double logH = newLL - thisTheta.LL;

	//accept or reject
	if(log(randomGenerator.getRand()) < logH){
		double* tmp = theseFrequencies;
		theseFrequencies = tmpBaseFreq;
		tmpBaseFreq = tmp;
		thisTheta.LL = newLL;
	}
};

void TThetaEstimatorRatio::oneMCMCIteration(TRandomGenerator & randomGenerator){
	//update theta1
	updateTheta(data, theta, baseFreq, theta2.logTheta - phiPriorMean, randomGenerator);

	//update theta2
	updateTheta(data2, theta2, baseFreq2, theta.logTheta + phiPriorMean, randomGenerator);

	//update base frequencies



};

