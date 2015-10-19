/*
 * TWindow.cpp
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#include "TWindow.h"

//-------------------------------------------------------
//EMConstants
//-------------------------------------------------------
EMParameters::EMParameters(){
	numIterations = -1;
	numThetaOnlyUpdates = -1;
	maxEpsilon = 0.0;
	NewtonRalphsonNumIterations = -1;
	NewtonRalphsonMaxF = 0.0;
	initalTheta = 0.0;
	initThetaSearchFactor = -1;
	initThetaNumSearchIterations = -1;
}

//-------------------------------------------------------
//Twindow
//-------------------------------------------------------
TWindow::TWindow(){
	start = -1;
	end = -1;
	length = -1;
	sites = NULL;
	sitesInitialized = false;
	coverage = -1.0;
	fractionSitesNoData = -1.0;
	fractionsitesCoverageAtLeastTwo = -1.0;
};

TWindow::TWindow(long Start, long End){
	start = Start;
	end = End;
	initSites(end - start); //end NOT in window!
};


void TWindow::clear(){
	for(int i=0; i<length; ++i) sites[i].clear();
	coverage = -1.0;
	fractionSitesNoData = -1.0;
	fractionsitesCoverageAtLeastTwo = -1.0;
};

void TWindow::move(long Start, long End){
	start = Start;
	end = End;
	if(sitesInitialized){
		if((start - end) != length)
			initSites(end - start);
		else  clear();
	} else initSites(end - start);
};

bool TWindow::addFromRead(BamTools::BamAlignment & bamAlignement, TReadGroups* readGroups){
	/* Note:
	 * Function returns true if read also maps to next window and
	 * returns false if end of read is within this (or a previous) window
	 */

	if(bamAlignement.Position >= end) return true;

	//find first position to be within window
	double len = bamAlignement.AlignedBases.length();
	if(bamAlignement.Position + len < start) return false;

	//find which position to consider first
	int firstPos = start - bamAlignement.Position;
	if(firstPos < 0) firstPos = 0;
	int lastPos = len;
	if(bamAlignement.Position + lastPos > end) lastPos = end - bamAlignement.Position;

	//Extract Read Group Info
	std::string readGroup;
	bamAlignement.GetTag("RG", readGroup);
	int readGroupId = readGroups->find(readGroup);

	//add sites
	int internalPos = bamAlignement.Position + firstPos - start;
	char base; BaseContext context;
	for(int pos = firstPos; pos < lastPos; ++pos, ++internalPos){
		/* Note:
		 * Reference is 5' -> 3'
		 * Hence for any read mapping "forward":
		 * 		- distance from 5' = pos + 1
		 * 		- distance from 3' = len - pos
		 * For any read mapping "reverse":
		 * 		- distance from 5' = len - pos
		 * 		- distance from 3' = pos + 1
		 * In both cases: 1) distance is 1-based!
		 *                2) Ignoring indels when calculating distances
		 *                3) Function add needs first pos5, then pos3
		 */

		//figure out context (base + previous base)
		base = bamAlignement.AlignedBases.at(pos);
		if(bamAlignement.IsReverseStrand()){
			if(pos >= bamAlignement.AlignedBases.length()) context = genoMap.getContext('N', base);
			else context = genoMap.getContext(bamAlignement.AlignedBases.at(pos + 1), base);

			sites[internalPos].add(base, bamAlignement.AlignedQualities.at(pos), pos, len - pos, pos + 1, context, readGroupId);
		} else {
			if(pos == 0) context = genoMap.getContext('N', base);
			else context = genoMap.getContext(bamAlignement.AlignedBases.at(pos - 1), base);
			sites[internalPos].add(base, bamAlignement.AlignedQualities.at(pos), pos, pos + 1, len - pos, context, readGroupId);
		}
	}

	//return if part of the read maps to next window
	if(lastPos == len) return false;
	else return true;
}

void TWindow::addReferenceBaseToSites(BamTools::Fasta & reference, int & refId){
	int stop = end - 1; //note that end is last position + 1
	std::string ref; //fasta object fills string
	reference.GetSequence(refId, start, stop-1, ref);

	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			sites[i].setRefBase(ref[i]);
		}
	}
}

void TWindow::applyMask(TBedReader* mask){
	//test if mask is required
	if(mask->hasPositionsInWindow(start)){
		//skip sites listed in mask by setting their hasData = false
		std::vector<long> thesePos = mask->getPositionInWindow(start);
		for(std::vector<long>::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
			sites[*it - start + 1].hasData = false;
		}
	}

	//also update counters!
}

void TWindow::estimateBaseFrequencies(){
	//estimate initial base frequencies
	baseFreq.clear();
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			sites[i].addToBaseFrequencies(baseFreq);
		}
	}
	baseFreq.normalize();
}

void TWindow::calculateEmissionProbabilities(TPMD & pmdObject, TRecalibration* recalObject){
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			recalObject->calcEmissionProbabilities(sites[i], pmdObject);
		}
	}
}

void TWindow::callMLEGenotype(TPMD & pmdObject, TRecalibration* recalObject, gz::ogzstream & out, std::string & chr, bool printAll){
	if(printAll){
		for(int i=0; i<length; ++i){
			out << chr << "\t" << start + i;
			if(sites[i].hasData) recalObject->calcEmissionProbabilities(sites[i], pmdObject);
			sites[i].callMLEGenotype(genoMap, out);
			out << "\n";
		}
	} else {
		for(int i=0; i<length; ++i){
			if(sites[i].hasData){
				out << chr << "\t" << start + i;
				recalObject->calcEmissionProbabilities(sites[i], pmdObject);
				sites[i].callMLEGenotype(genoMap, out);
				out << "\n";
			}
		}
	}
}

void TWindow::printPileup(TPMD & pmdObject, TRecalibration* recalObject, std::ofstream & out, std::string & chr){
	//calc emission probs
	for(int i=0; i<length; ++i){
		recalObject->calcEmissionProbabilities(sites[i], pmdObject);
	}
	//print pileup
	for(int i=0; i<length; ++i){
		out << chr << "\t" << start + i + 1 << "\t" << sites[i].bases.size() << "\t" << sites[i].getBases() << "\t" << sites[i].getEmissionProbs() << "\n";
	}
}

void TWindow::calcCoverage(){
	//calculate and return coverage
	coverage = 0.0;
	long noData = 0;
	long plentyData = 0;
	for(int i=0; i<length; ++i){
		coverage += sites[i].bases.size();
		if(sites[i].bases.size() == 0) ++ noData;
		else if(sites[i].bases.size() > 1) ++ plentyData;
	}

	coverage = coverage / (double) length;
	fractionSitesNoData = (double) noData / (double) length;
	fractionsitesCoverageAtLeastTwo = (double) plentyData / (double) length;
}

double TWindow::calcLogLikelihood(double* pGenotype){
	double LL = 0.0;
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			LL += sites[i].calculateLogLikelihood(pGenotype);
		}
	}
	return LL;
}

void TWindow::addSitesToBQSR(TRecalibrationBQSR & bqsr, TLog* logfile){
	logfile->listFlush("Adding sites to BQSR ...");
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			bqsr.addSite(sites[i]);
		}
	}
	logfile->write(" done!");
}
//-------------------------------------------------------
//TwindowDiploid
//-------------------------------------------------------
void TWindowDiploid::initSites(long newLength){
	if(sitesInitialized)
		delete[] sites;
	length = newLength;
	sites = new TSiteDiploid[length];
	sitesInitialized = true;
	coverage = -1.0;
	fractionSitesNoData = -1.0;
	fractionsitesCoverageAtLeastTwo = -1.0;
}

void TWindowDiploid::fillPGenotype(double* pGenotype, double & expTheta){
	for(int i=0; i<4; ++i){
		//homozygous genotypes
		pGenotype[genoMap.getGenotype(i,i)] = baseFreq[i] * (expTheta + baseFreq[i] * (1.0 - expTheta));
		//heterozygous genotypes
		for(int j=i+1; j<4; ++j){
			pGenotype[genoMap.getGenotype(i,j)] = 2.0 * baseFreq[i] * baseFreq[j] *  (1.0 - expTheta);
		}
	}
}

void TWindowDiploid::fillP_G(double* P_G, double* pGenotype){
	for(int g=0; g<10; ++g)
		P_G[g] = 0.0;

	//calculate P_g for each site
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			sites[i].calculateP_g(pGenotype);
			for(int g=0; g<10; ++g){
				P_G[g] += sites[i].P_g[g];
			}
		}
	}
}

void TWindowDiploid::estimateTheta(EMParameters & EMParams, TPMD & pmdObject, TRecalibration* recalObject, std::ofstream & out, TLog* logfile){
	logfile->startIndent("Estimating Theta:");

	//measure runtime
	struct timeval startTime, endTime;
	gettimeofday(&startTime, NULL);

	//Theta container
	Theta thetaContainer;

	//estimate initial base frequencies
	//calculate per site emission probabilities
	logfile->listFlush("Calculating emission probabilities ...");
	calculateEmissionProbabilities(pmdObject, recalObject);
	logfile->write(" done!");

	//get num sites with data
	int lengthWithData = 0;
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			++lengthWithData;
		}
	}

	//estimate starting parameters
	logfile->startIndent("Estimating initial parameters:");
	logfile->listFlush("Estimating initial base frequencies ...");
	estimateBaseFrequencies();
	logfile->write(" done!");
	logfile->conclude("Pi(A) = " + toString(baseFreq[0]) + ", Pi(C) = " + toString(baseFreq[1]) + ", Pi(G) = " + toString(baseFreq[2]) + ", Pi(T) = " + toString(baseFreq[3]));

	//set initial parameters
	logfile->listFlush("Estimating initial theta ...");
	findGoodStartingTheta(thetaContainer, EMParams);
	logfile->write(" done!");
	logfile->conclude("Starting EM with theta = ", thetaContainer.theta);
	logfile->endIndent();

	//Run EM
	logfile->listFlush("Running EM to find ML estimations ...");
	runEMForTheta(thetaContainer, EMParams, lengthWithData);
	logfile->write(" done!");
	logfile->conclude("theta was estimated at ", thetaContainer.theta);

	//confidence intervals
	logfile->listFlush("Estimating approximate confidence intervals from Fisher-Information ...");
	estimateConfidenceInterval(thetaContainer);
	logfile->write(" done!");
	logfile->conclude("95% confidence intervals are theta +- " + toString(thetaContainer.thetaConfidence));

	//write results to file
	//position
	out << start << "\t" << end-1;
	//coverage NOTE: assumes coverage has been calculated before...
	out << "\t" << coverage << "\t" << fractionSitesNoData << "\t" << fractionsitesCoverageAtLeastTwo;
	//estimated params
	for(int i=0; i<4; ++i)
		out << "\t" << baseFreq[i];
	out << "\t" << thetaContainer.theta << "\t" << thetaContainer.theta - thetaContainer.thetaConfidence << "\t" << thetaContainer.theta + thetaContainer.thetaConfidence << "\t" << thetaContainer.LL << std::endl;

	//finish
	gettimeofday(&endTime, NULL);
	logfile->list("Total computation time for this window was ", endTime.tv_sec  - startTime.tv_sec, "s");
	logfile->endIndent();
}


void TWindowDiploid::findGoodStartingTheta(Theta & thetaContainer, EMParameters & EMParams){
	//assumes that initial base frequencies have been estimated and site emission probs calculated!
	double pGenotype[10];
	double initTheta = EMParams.initalTheta;
	double oldTheta = initTheta;
	double expTheta = exp(-initTheta);

	//calc initial LL
	fillPGenotype(pGenotype, expTheta);
	thetaContainer.LL = calcLogLikelihood(pGenotype);

	//run iterations
	double oldLL = thetaContainer.LL;
	double factor = EMParams.initThetaSearchFactor;
	int numUpdates;
	for(int i=0; i<EMParams.initThetaNumSearchIterations; ++i){
		//first test increase in theta
		numUpdates = -1;
		do{
			++numUpdates;
			oldLL = thetaContainer.LL;
			oldTheta = initTheta;
			initTheta *= factor;
			expTheta = exp(-initTheta);
			fillPGenotype(pGenotype, expTheta);
			thetaContainer.LL = calcLogLikelihood(pGenotype);
		} while(oldLL < thetaContainer.LL);
		if(numUpdates == 0){
			//then test decrease in theta
			initTheta = oldTheta;
			thetaContainer.LL = oldLL;
			//maybe smaller?
			do{
				oldLL = thetaContainer.LL;
				oldTheta = initTheta;
				initTheta /= factor;
				expTheta = exp(-initTheta);
				fillPGenotype(pGenotype, expTheta);
				thetaContainer.LL = calcLogLikelihood(pGenotype);
			} while(oldLL < thetaContainer.LL);
		}
		factor = sqrt(factor);
		initTheta = oldTheta;
		thetaContainer.LL = oldLL;
	}
	//return previous
	thetaContainer.setTheta(oldTheta);
	thetaContainer.LL = oldLL;

	//check if values make sense. If theta < 1/(10*windowsize), set it to 1/(10*windowsize)
	if(thetaContainer.theta < 0.1/length){
		thetaContainer.setTheta(0.1/length);
	} else if(thetaContainer.theta > 1.0){
		thetaContainer.setTheta(1.0);
	}
}

void TWindowDiploid::runEMForTheta(Theta & thetaContainer, EMParameters & EMParams, int & lengthWithData){
	//prepare storage
	double pGenotype[10];
	double P_G[10];
	double tmp[4];
	double tmpSum;
	arma::mat Jacobian(6,6);
	arma::vec F(6);
	arma::mat JxF(6,6);
	Genotype geno;
	double maxF;
	int failedAttempts = 0;
	double oldTheta, rho, mu = lengthWithData;
	double oldLL = -9e100;

	//start EM loop
	int numThetaOnlyUpdatesDone = EMParams.numThetaOnlyUpdates; //do regular step first
	numThetaOnlyUpdatesDone = 0;
	int totIterations = EMParams.numIterations * EMParams.numThetaOnlyUpdates;
	for(int iter = 0; iter < totIterations; ++iter){
		//a) pre-calc expTheta
		oldTheta = thetaContainer.theta;
		rho = thetaContainer.expTheta / (1.0 - thetaContainer.expTheta);

		//b) calculate	substitution probabilities
		fillPGenotype(pGenotype, thetaContainer.expTheta);

		//c) Calculate all genotype probabilities for all sites
		fillP_G(P_G, pGenotype);

		//d) Find new parameter estimates using Newton-Ralphson
		if(numThetaOnlyUpdatesDone < EMParams.numThetaOnlyUpdates){
			//update only theta: most difficult parameter and it is much faster to update only this one alone.
			for(int n=0; n<EMParams.NewtonRalphsonNumIterations; ++n){
				//i) calculate F() (Note: index is zero based!)
				F(4) = lengthWithData;
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
				if(F(4) < EMParams.NewtonRalphsonMaxF){
					thetaContainer.setTheta(-log(rho / (1.0 + rho)));
					break;
				}
			}
			++numThetaOnlyUpdatesDone;
			if(thetaContainer.theta == oldTheta) numThetaOnlyUpdatesDone = EMParams.numThetaOnlyUpdates;
		} else {
			numThetaOnlyUpdatesDone = 0;
			//update all parameters in EM
			for(int n=0; n<EMParams.NewtonRalphsonNumIterations; ++n){
				//i) calculate F (Note: index is zero based!)
				F(4) = lengthWithData;
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
					if(maxF < EMParams.NewtonRalphsonMaxF){
						thetaContainer.setTheta(-log(rho / (1.0 + rho)));
						break;
					}
				} else {
					++failedAttempts;
					//solve did not work -> start with higher theta!
					thetaContainer.setTheta(EMParams.initalTheta);
					for(int i=0; i<failedAttempts; ++i)
						thetaContainer.theta *= 10.0;
					//reset others
					mu = lengthWithData;
					thetaContainer.LL = -9e100;
					iter = 0;
					numThetaOnlyUpdatesDone = 0;
					break;
				}
			}
		}

		//e) do we break EM? Check LL
		if(iter > 0 && iter % EMParams.numThetaOnlyUpdates == 0){
			oldLL = thetaContainer.LL;
			thetaContainer.LL = calcLogLikelihood(pGenotype);
			if(thetaContainer.LL > -9e100 && (thetaContainer. LL - oldLL) < EMParams.maxEpsilon) break;

			//maybe theta = 0?
			if(thetaContainer.theta < 0.1/length){
				oldLL = thetaContainer.LL;
				oldTheta = thetaContainer.theta;
				//test with theta = 0.0
				thetaContainer.setTheta(0.0);
				fillPGenotype(pGenotype, thetaContainer.expTheta);
				thetaContainer.LL = calcLogLikelihood(pGenotype);

				if(thetaContainer.LL < oldLL){
					thetaContainer.setTheta(oldTheta);
					thetaContainer.LL = oldLL;
				}
				break;
			}
		}

		//For debugging
		//std::cout << std::setprecision(9) << iter << ") theta = " << theta << "\tLL = " << LL << "\teps = " << fabs(oldLL - LL) << std::endl;
	}
}

void TWindowDiploid::estimateConfidenceInterval(Theta & thetaContainer){
	//we estimate an approximate confidence interval for theta using the Fisher information
	//This function assumes that EM has already been run!

	//calculate P(g|theta, pi)
	double pGenotype[10];
	fillPGenotype(pGenotype, thetaContainer.expTheta);

	//calclate d/dtheta P(g|theta, pi)
	double deriv_pGenotype[10];
	for(int k=0; k<4; ++k){
		//homozygous genotype
		deriv_pGenotype[genoMap.getGenotype(k, k)] = (baseFreq[k] * baseFreq[k] - baseFreq[k]) * thetaContainer.expTheta;
		//heterozygous genotypes
		for(int l=k+1; l<4; ++l){
			deriv_pGenotype[genoMap.getGenotype(k, l)] = 2.0 * baseFreq[k] * baseFreq[l] * thetaContainer.expTheta;
		}
	}

	//sum Ri over all sites
	double FisherInfo = 0.0;
	double Ri;
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			//calc Ri
			Ri = sites[i].calculateWeightedSumOfEmissionProbs(deriv_pGenotype) / sites[i].calculateWeightedSumOfEmissionProbs(pGenotype);
			//add to Fisher Info
			FisherInfo += Ri * (Ri + 1.0);
		}
	}

	//estimate confidence interval
	thetaContainer.thetaConfidence = 1.96 / sqrt(FisherInfo);
}



void TWindowDiploid::calcLikelihoodSurface(TPMD & pmdObject, TRecalibration* recalObject, std::ofstream & out, int & steps){
	//estimate initial base frequencies
	//calculate per site emission probabilities
	calculateEmissionProbabilities(pmdObject, recalObject);
	estimateBaseFrequencies();

	//write header
	out << "log10(theta)\ttheta\tLL\n";

	//prepare storage
	double pGenotype[10];

	//calculate likelihood surface
	double minLogTheta = -5.0;
	double maxLogTheta = 2.0;
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
		fillPGenotype(pGenotype, expTheta);
		LL = calcLogLikelihood(pGenotype);

		//write results
		out << logTheta << "\t" << theta << "\t" << LL << "\n";
	}
}


//-------------------------------------------------------
//TWindowHaploid
//-------------------------------------------------------

void TWindowHaploid::initSites(long newLength){
	if(sitesInitialized)
		delete[] sites;
	length = newLength;
	sites = new TSiteHaploid[length];
	sitesInitialized = true;
	coverage = -1.0;
	fractionSitesNoData = -1.0;
	fractionsitesCoverageAtLeastTwo = -1.0;
}

void TWindowHaploid::fillPGenotype(double* pGenotype){
	for(int i=0; i<4; ++i){
		pGenotype[i] = baseFreq[i];
	}
}

double TWindowHaploid::calcLogLikelihood(){
	double pGenotype[4];
	fillPGenotype(pGenotype);

	double LL = 0.0;
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			LL += sites[i].calculateLogLikelihood(pGenotype);
		}
	}
	return LL;
}

void TWindowHaploid::addToJacobian(TRecalibrationEM* reclObject){
	//assumes that frequencies have been calculated!
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			reclObject->addSiteToJacobianAndF(sites[i].bases, &baseFreq);
		}
	}
}

void TWindowHaploid::addToLikelihoodRecalibration(TRecalibrationEM* reclObject){
	//assumes that frequencies have been calculated!
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			reclObject->addSiteToLikelihood(sites[i].bases, &baseFreq);
		}
	}
}
