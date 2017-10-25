/*
 * TThetaEstimator.cpp
 *
 *  Created on: Sep 25, 2017
 *      Author: phaentu
 */

#include "TThetaEstimator.h"

//-------------------------------------------------------
//TThetaEstimator
//-------------------------------------------------------
TThetaEstimator::TThetaEstimator(TParameters & params, TLog* Logfile){
	logfile = Logfile;
	numGenotypes = 10;

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

	//counters and tmp variables
	init();

	//done
	logfile->endIndent();
}

TThetaEstimator::TThetaEstimator(TLog* Logfile){
	logfile = Logfile;
	numGenotypes = 10;

	//set EM parameters to default
	numIterations = -1;
	numThetaOnlyUpdates = -1;
	maxEpsilon = 0.0;
	NewtonRaphsonNumIterations = -1;
	NewtonRaphsonMaxF = 0.0;
	initalTheta = 0.0;
	initThetaSearchFactor = -1;
	initThetaNumSearchIterations = -1;

	init();
}

void TThetaEstimator::init(){
	//initialize arrays
	pGenotype = new double[10];
	P_G = new double[10];
	baseFreq = new double[4];

	//empty
	clear();

	//tmp stuff
	g = 0;
	doublePointer = NULL;
	sum = 0.0;
}

void TThetaEstimator::add(TSite & site){
	//assumes that emission probabilities were calculated!!
	++totNumSitesAdded;
	cumulativeDepth += site.depth();

	//add if site has data
	if(site.hasData){
		//store emission probabilities
		sites.push_back(new double[10]);
		doublePointer = *sites.rbegin();

		for(g=0; g<numGenotypes; ++g)
			doublePointer[g] = site.emissionProbabilities[g];

		++numSitesWithData;

		//add site to base frequency estimation
		site.addToBaseFrequencies(initialBaseFreq);

		//count sites covered >=2
		if(site.depth() > 1)
			++numSitesCoveredTwiceOrMore;
	}
}

void TThetaEstimator::fillPGenotype(double* & pGeno, double & expTheta){
	//assumes that base frequencies are set!
	for(int i=0; i<4; ++i){
		//homozygous genotypes
		pGeno[genoMap.getGenotype(i,i)] = baseFreq[i] * (expTheta + baseFreq[i] * (1.0 - expTheta));
		//heterozygous genotypes
		for(int j=i+1; j<4; ++j){
			pGeno[genoMap.getGenotype(i,j)] = 2.0 * baseFreq[i] * baseFreq[j] *  (1.0 - expTheta);
		}
	}
}

void TThetaEstimator::fillPGenotype(double & expTheta){
	fillPGenotype(pGenotype, expTheta);
}

void TThetaEstimator::fillPGenotype(double* & pGeno){
	double expTheta = exp(-theta.theta);
	fillPGenotype(pGeno, expTheta);
}

void TThetaEstimator::fillP_G(std::vector<double*> & theseSites){
	//assumes that pGenotype is set!
	for(g=0; g<numGenotypes; ++g)
		P_G[g] = 0.0;

	//calculate P_g for each site
	for(siteIt=theseSites.begin(); siteIt != theseSites.end(); ++siteIt){
		sum = 0.0;
		for(g=0; g<numGenotypes; ++g){
			P_g_oneSite[g] =  (*siteIt)[g] * pGenotype[g];
			sum += P_g_oneSite[g];
		}
		for(g=0; g<numGenotypes; ++g)
			P_G[g] += P_g_oneSite[g] / sum;
	}
}

double TThetaEstimator::calcLogLikelihood(std::vector<double*> & theseSites){
	double LL = 0.0;
	for(siteIt=sites.begin(); siteIt != sites.end(); ++siteIt){
		sum = 0.0;
		for(g=0; g<numGenotypes; ++g)
			sum +=  (*siteIt)[g] * pGenotype[g];
		LL += log(sum);
	}
	return LL;
}

void TThetaEstimator::findGoodStartingTheta(std::vector<double*> & theseSites){
	//set base frequencies to initial base frequencies
	initialBaseFreq.normalize();
	for(int i=0; i<4; ++i)
		baseFreq[i] = initialBaseFreq[i];

	//variables
	double initTheta = initalTheta;
	double oldTheta = initTheta;
	double expTheta = exp(-initTheta);

	//calc initial LL
	fillPGenotype(expTheta);
	theta.LL = calcLogLikelihood(theseSites);

	//run iterations
	double oldLL = theta.LL;
	double factor = initThetaSearchFactor;
	int numUpdates;
	for(int i=0; i<initThetaNumSearchIterations; ++i){
		//first test increase in theta
		numUpdates = -1;
		do{
			++numUpdates;
			oldLL = theta.LL;
			oldTheta = initTheta;
			initTheta *= factor;
			expTheta = exp(-initTheta);
			fillPGenotype(expTheta);
			theta.LL = calcLogLikelihood(theseSites);
		} while(oldLL < theta.LL);
		if(numUpdates == 0){
			//then test decrease in theta
			initTheta = oldTheta;
			theta.LL = oldLL;
			//maybe smaller?
			do{
				oldLL = theta.LL;
				oldTheta = initTheta;
				initTheta /= factor;
				expTheta = exp(-initTheta);
				fillPGenotype(expTheta);
				theta.LL = calcLogLikelihood(theseSites);
			} while(oldLL < theta.LL);
		}
		factor = sqrt(factor);
		initTheta = oldTheta;
		theta.LL = oldLL;
	}
	//return previous
	theta.setTheta(oldTheta);
	theta.LL = oldLL;

	//check if values make sense. If theta < 1/(10*windowsize), set it to 1/(10*windowsize)
	if(theta.theta < 0.1/(double) totNumSitesAdded){
		theta.setTheta(0.1/(double) totNumSitesAdded);
	} else if(theta.theta > 1.0){
		theta.setTheta(1.0);
	}
}

void TThetaEstimator::runEMForTheta(std::vector<double*> & theseSites){
	//prepare storage
	double tmp[4];
	double tmpSum;
	arma::mat Jacobian(6,6);
	arma::vec F(6);
	arma::mat JxF(6,6);
	Genotype geno;
	double maxF;
	int failedAttempts = 0;
	double oldTheta, rho, mu = numSitesWithData;
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
		fillPGenotype(theta.expTheta);

		//c) Calculate all genotype probabilities for all sites
		fillP_G(theseSites);

		//d) Find new parameter estimates using Newton-Raphson
		if(numThetaOnlyUpdatesDone < numThetaOnlyUpdates){
			//update only theta: most difficult parameter and it is much faster to update only this one alone.
			for(int n=0; n<NewtonRaphsonNumIterations; ++n){
				//i) calculate F() (Note: index is zero based!)
				F(4) = numSitesWithData;
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
				F(4) = numSitesWithData;
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
					mu = numSitesWithData;
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
			theta.LL = calcLogLikelihood(theseSites);
			if(theta.LL > -9e100 && (theta.LL - oldLL) < maxEpsilon) break;

			//maybe theta = 0?
			if(theta.theta < 0.1/(double) totNumSitesAdded){
				oldLL = theta.LL;
				oldTheta = theta.theta;

				//test with theta = 0.0
				theta.setTheta(0.0);
				fillPGenotype(theta.expTheta);
				theta.LL = calcLogLikelihood(theseSites);

				if(theta.LL < oldLL){
					theta.setTheta(oldTheta);
					theta.LL = oldLL;
				}
				break;
			}
		}

		//For debugging
		//std::cout << std::setprecision(9) << iter << ") theta = " << thetaContainer.theta << "\tLL = " << thetaContainer.LL << "\teps = " << fabs(oldLL - thetaContainer.LL) << std::endl;
	}
}

void TThetaEstimator::estimateConfidenceInterval(std::vector<double*> & theseSites){
	//we estimate an approximate confidence interval for theta using the Fisher information
	//This function assumes that EM has already been run!

	//calculate P(g|theta, pi)
	fillPGenotype(theta.expTheta);

	//calclate d/dtheta P(g|theta, pi)
	double deriv_pGenotype[10];
	for(int k=0; k<4; ++k){
		//homozygous genotype
		deriv_pGenotype[genoMap.getGenotype(k, k)] = (baseFreq[k] * baseFreq[k] - baseFreq[k]) * theta.expTheta;
		//heterozygous genotypes
		for(int l=k+1; l<4; ++l){
			deriv_pGenotype[genoMap.getGenotype(k, l)] = 2.0 * baseFreq[k] * baseFreq[l] * theta.expTheta;
		}
	}

	//sum Ri over all sites
	double FisherInfo = 0.0;
	double Ri, Ri_a, Ri_b;

	for(siteIt=theseSites.begin(); siteIt != theseSites.end(); ++siteIt){
		//calc Ri
		Ri_a = 0.0; Ri_b = 0.0;
		for(g=0; g<numGenotypes; ++g){
			Ri_a += (*siteIt)[g] * deriv_pGenotype[g];
			Ri_b += (*siteIt)[g] * pGenotype[g];
		}
		Ri = Ri_a / Ri_b;

		//add to Fisher Info
		FisherInfo += Ri * (Ri + 1.0);
	}

	//estimate confidence interval
	//TODO: Fisher Info can be negative -> SQRT will be nan!
	theta.thetaConfidence = 1.96 / sqrt(FisherInfo);
}

bool TThetaEstimator::estimateTheta(std::vector<double*> & theseSites){
	if(numSitesWithData < 100){
		logfile->write("Can not estimate theta, less than 100 sites with data in this region!");
		return false;
	}

	//estimate starting parameters
	logfile->listFlush("Estimating initial parameters ...");
	findGoodStartingTheta(theseSites);
	logfile->done();
	logfile->conclude("Initial base frequencies: Pi(A) = " + toString(baseFreq[0]) + ", Pi(C) = " + toString(baseFreq[1]) + ", Pi(G) = " + toString(baseFreq[2]) + ", Pi(T) = " + toString(baseFreq[3]));
	logfile->conclude("Starting EM with theta = ", theta.theta);

	//Run EM
	logfile->listFlush("Running EM to find ML estimate ...");
	runEMForTheta(theseSites);
	logfile->done();
	logfile->conclude("theta was estimated at ", theta.theta);

	//confidence intervals
	logfile->listFlush("Estimating approximate confidence intervals from Fisher-Information ...");
	estimateConfidenceInterval(theseSites);
	logfile->done();
	logfile->conclude("95% confidence intervals are theta +- " + toString(theta.thetaConfidence));
	return true;
}

//------------------------------------------------------------
//Functions to run estimation-
//------------------------------------------------------------
bool TThetaEstimator::estimateTheta(){
	return estimateTheta(sites);
}

void TThetaEstimator::setTheta(double Theta){
	theta.setTheta(Theta);
}

void TThetaEstimator::setBaseFreq(TBaseFrequencies & BaseFreq){
	for(int i=0; i<4; ++i)
		baseFreq[i] = BaseFreq[i];
}

void TThetaEstimator::writeHeader(std::ofstream & out){
	out << "depth\tfracMissing\tfracTwoOrMore\tpi(A)\tpi(C)\tpi(G)\tpi(T)\ttheta_MLE\ttheta_C95_l\ttheta_C95_u\tLL";
}

void TThetaEstimator::writeResultsToFile(std::ofstream & out){
	out << "\t" << cumulativeDepth / (double) totNumSitesAdded << "\t" << (double) (totNumSitesAdded - numSitesWithData) / (double) totNumSitesAdded << "\t" << (double) numSitesCoveredTwiceOrMore / (double) totNumSitesAdded;	//estimated params
	for(int i=0; i<4; ++i)
		out << "\t" << baseFreq[i];

	out << "\t" << theta.theta;
	out << "\t" << theta.theta - theta.thetaConfidence;
	out << "\t" << theta.theta + theta.thetaConfidence;
	out << "\t" << theta.LL;
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
		fillPGenotype(expTheta);
		LL = calcLogLikelihood(sites);

		//write results
		out << std::setprecision(12) << logTheta << "\t" << theta << "\t" << LL << "\n";
	}
}

void TThetaEstimator::bootstrapTheta(TRandomGenerator & randomGenerator, std::ofstream & out){
	logfile->listFlush("Bootstrapping sites ...");

	//prepare variables
	std::vector<double*> bootstrappedSites;
	int r = 0;

	//how many sites with data will we have? Draw from binomial
	double probHasData = (double) numSitesWithData / (double) totNumSitesAdded;
	long numBootstrappedSites = randomGenerator.getBiomialRand(probHasData, sites.size());

	//now pick among sites with data with replacment
	for(long i=0; i<numBootstrappedSites; ++i){
		r = 0;
		while(r == 0) r = randomGenerator.pickOne(numBootstrappedSites - 1);
		bootstrappedSites.push_back(sites.at(r));
	}
	logfile->done();

	//estimate theta
	estimateTheta(bootstrappedSites);

	//write output (modified from standard output
	out << "\t" << (double) (totNumSitesAdded - numBootstrappedSites) / (double) totNumSitesAdded;	//estimated params
	for(int i=0; i<4; ++i)
		out << "\t" << baseFreq[i];
	out << "\t" << theta.theta << "\t" << theta.theta - theta.thetaConfidence << "\t" << theta.theta + theta.thetaConfidence << "\t" << theta.LL;

	//clean up
	bootstrappedSites.clear();
}


