/*
 * TThetaEstimator.cpp
 *
 *  Created on: Sep 25, 2017
 *      Author: phaentu
 */

#include "TThetaEstimator.h"

void TThetaEstimator::add(TSite & site){
	//assumes that emission probabilities were calculated!!
	++totNumSitesAdded;
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

void TThetaEstimator::fillPGenotype(double & expTheta){
	//assumes that base frequencies are set!
	for(int i=0; i<4; ++i){
		//homozygous genotypes
		pGenotype[genoMap.getGenotype(i,i)] = baseFreq[i] * (expTheta + baseFreq[i] * (1.0 - expTheta));
		//heterozygous genotypes
		for(int j=i+1; j<4; ++j){
			pGenotype[genoMap.getGenotype(i,j)] = 2.0 * baseFreq[i] * baseFreq[j] *  (1.0 - expTheta);
		}
	}
}

void TThetaEstimator::fillP_G(){
	//assumes that pGenotype is set!
	for(g=0; g<numGenotypes; ++g)
		P_G[g] = 0.0;

	//calculate P_g for each site
	for(siteIt=sites.begin(); siteIt != sites.end(); ++siteIt){
		sum = 0.0;
		for(g=0; g<numGenotypes; ++g){
			P_g_oneSite[g] =  (*siteIt)[g] * pGenotype[g];
			sum += P_g_oneSite[g];
		}
		for(g=0; g<numGenotypes; ++g)
			P_G[g] += P_g_oneSite[g] / sum;
	}
}

double TThetaEstimator::calcLogLikelihood(){
	double LL = 0.0;
	for(siteIt=sites.begin(); siteIt != sites.end(); ++siteIt){
		sum = 0.0;
		for(g=0; g<numGenotypes; ++g)
			sum +=  (*siteIt)[g] * pGenotype[g];
		LL += log(sum);
	}
	return LL;
}

void TThetaEstimator::findGoodStartingTheta(EMParameters & EMParams){
	//assumes that initial base frequencies have been estimated
	double initTheta = EMParams.initalTheta;
	double oldTheta = initTheta;
	double expTheta = exp(-initTheta);

	//calc initial LL
	fillPGenotype(expTheta);
	theta.LL = calcLogLikelihood();

	//run iterations
	double oldLL = theta.LL;
	double factor = EMParams.initThetaSearchFactor;
	int numUpdates;
	for(int i=0; i<EMParams.initThetaNumSearchIterations; ++i){
		//first test increase in theta
		numUpdates = -1;
		do{
			++numUpdates;
			oldLL = theta.LL;
			oldTheta = initTheta;
			initTheta *= factor;
			expTheta = exp(-initTheta);
			fillPGenotype(expTheta);
			theta.LL = calcLogLikelihood();
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
				theta.LL = calcLogLikelihood();
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


void TThetaEstimator::runEMForTheta(EMParameters & EMParams){
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
	int numThetaOnlyUpdatesDone = EMParams.numThetaOnlyUpdates; //do regular step first
	numThetaOnlyUpdatesDone = 0;
	int totIterations = EMParams.numIterations * EMParams.numThetaOnlyUpdates;
	for(int iter = 0; iter < totIterations; ++iter){
		//a) pre-calc expTheta
		oldTheta = theta.theta;
		rho = theta.expTheta / (1.0 - theta.expTheta);

		//b) calculate	substitution probabilities
		fillPGenotype(theta.expTheta);

		//c) Calculate all genotype probabilities for all sites
		fillP_G();

		//d) Find new parameter estimates using Newton-Raphson
		if(numThetaOnlyUpdatesDone < EMParams.numThetaOnlyUpdates){
			//update only theta: most difficult parameter and it is much faster to update only this one alone.
			for(int n=0; n<EMParams.NewtonRaphsonNumIterations; ++n){
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
				if(F(4) < EMParams.NewtonRaphsonMaxF){
					theta.setTheta(-log(rho / (1.0 + rho)));
					break;
				}
			}
			++numThetaOnlyUpdatesDone;
			if(theta.theta == oldTheta) numThetaOnlyUpdatesDone = EMParams.numThetaOnlyUpdates;
		} else {
			numThetaOnlyUpdatesDone = 0;
			//update all parameters in EM
			for(int n=0; n<EMParams.NewtonRaphsonNumIterations; ++n){
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

					if(maxF < EMParams.NewtonRaphsonMaxF || n == (EMParams.NewtonRaphsonNumIterations-1)){
						theta.setTheta(-log(rho / (1.0 + rho)));
						break;
					}
				} else {
					++failedAttempts;

					//solve did not work -> start with higher theta!
					theta.setTheta(EMParams.initalTheta);
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
		if(iter > 0 && iter % EMParams.numThetaOnlyUpdates == 0){
			oldLL = theta.LL;
			theta.LL = calcLogLikelihood();
			if(theta.LL > -9e100 && (theta.LL - oldLL) < EMParams.maxEpsilon) break;

			//maybe theta = 0?
			if(theta.theta < 0.1/(double) totNumSitesAdded){
				oldLL = theta.LL;
				oldTheta = theta.theta;
				//test with theta = 0.0
				theta.setTheta(0.0);
				fillPGenotype(theta.expTheta);
				theta.LL = calcLogLikelihood();

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

void TThetaEstimator::estimateConfidenceInterval(){
	//we estimate an approximate confidence interval for theta using the Fisher information
	//This function assumes that EM has already been run!

	//calculate P(g|theta, pi)
	double pGenotype[10];
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
	for(siteIt=sites.begin(); siteIt != sites.end(); ++siteIt){
		//calc Ri
		Ri_a = 0.0; Ri_b = 0.0;
		for(g=0; g<numGenotypes; ++g){
			Ri_a += *siteIt[g] * deriv_pGenotype[g];
			Ri_a += *siteIt[g] * pGenotype[g];
		}
		Ri = Ri_a / Ri_b;

		//add to Fisher Info
		FisherInfo += Ri * (Ri + 1.0);
	}

	//estimate confidence interval
	theta.thetaConfidence = 1.96 / sqrt(FisherInfo);
}

//------------------------------------------------------------
//Public Functons
//------------------------------------------------------------

void TThetaEstimator::estimateTheta(EMParameters & EMParams, TRecalibration* recalObject, std::ofstream & out, TLog* logfile, bool & considerRegions){
	//estimate starting parameters
	logfile->startIndent("Estimating initial parameters:");
	logfile->listFlush("Initial base frequencies: Pi(A) = " + toString(baseFreq[0]) + ", Pi(C) = " + toString(baseFreq[1]) + ", Pi(G) = " + toString(baseFreq[2]) + ", Pi(T) = " + toString(baseFreq[3]));

	//set initial parameters
	logfile->listFlush("Estimating initial theta ...");
	findGoodStartingTheta(EMParams);
	logfile->write(" done!");
	logfile->conclude("Starting EM with theta = ", theta.theta);
	logfile->endIndent();

	//Run EM
	logfile->listFlush("Running EM to find ML estimate ...");
	runEMForTheta(EMParams);
	logfile->write(" done!");
	logfile->conclude("theta was estimated at ", thetaContainer.theta);

	//confidence intervals
	logfile->listFlush("Estimating approximate confidence intervals from Fisher-Information ...");
	estimateConfidenceInterval(thetaContainer);
	logfile->write(" done!");
	logfile->conclude("95% confidence intervals are theta +- " + toString(thetaContainer.thetaConfidence));

	//write results to file
	//position
	//TODO: think about how to treat when regions are to be considered
	if(considerRegions) out << "givenByUser\tgivenByUser\tgivenByUser";
	else out << start << "\t" << end-1;
	//coverage NOTE: assumes coverage has been calculated before...
	if(considerRegions) out << "\t" << "NA" << "\t" << "NA" << "\t" << "NA";
	else out << "\t" << coverage << "\t" << fractionSitesNoData << "\t" << fractionsitesCoverageAtLeastTwo;
	//estimated params
	for(int i=0; i<4; ++i)
		out << "\t" << baseFreq[i];
	out << "\t" << thetaContainer.theta << "\t" << thetaContainer.theta - thetaContainer.thetaConfidence << "\t" << thetaContainer.theta + thetaContainer.thetaConfidence << "\t" << thetaContainer.LL << std::endl;

	//finish
	gettimeofday(&endTime, NULL);
	logfile->list("Total computation time for this window was ", endTime.tv_sec  - startTime.tv_sec, "s");
	logfile->endIndent();
}

void TThetaEstimator::calcLikelihoodSurface(std::ofstream & out, int & steps){
	//write header
	out << "log10(theta)\ttheta\tLL\n";

	//prepare storage
	double pGenotype[10];

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
		LL = calcLogLikelihood();

		//write results
		out << std::setprecision(12) << logTheta << "\t" << theta << "\t" << LL << "\n";
	}
}



