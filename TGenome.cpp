/*
 * TGenome.cpp
 *
 *  Created on: Mar 15, 2015
 *      Author: wegmannd
 */


#include "TGenome.h"

//-------------------------------------------------------
//Constants
//-------------------------------------------------------
void Constants::initPmd(TParameters & params, TLog* logfile){
	if(params.parameterExists("pmd")){
		std::string pmd = params.getParameterString("pmd");
		logfile->list("Initializing post mortem damage for both C->T and G->A with function '" + pmd +"'.");
		pmd_C_T.initializeFunction(pmd);
		pmd_G_A.initializeFunction(pmd);
		logfile->conclude(pmd_C_T.getFunctionString());
	} else {
		//separate model for C->T and G->A
		//first C->T
		if(!params.parameterExists("pmdCT")) throw "Problem initializing post mortem damage: argument 'pmd' or 'pmdCT' has to be provided!";
		std::string pmd = params.getParameterString("pmdCT");
		logfile->list("Initializing post mortem C->T damage with function '" + pmd +"'.");
		pmd_C_T.initializeFunction(pmd);
		logfile->conclude(pmd_C_T.getFunctionString());

		//second G->A
		if(!params.parameterExists("pmdGA")) throw "Problem initializing post mortem damaga: argument 'pmd' or 'pmdGA' has to be provided!";
		pmd = params.getParameterString("pmdGA");
		logfile->list("Initializing post mortem G->A damage with function '" + pmd +"'.");
		pmd_G_A.initializeFunction(pmd);
		logfile->conclude(pmd_G_A.getFunctionString());
	}
}

//-------------------------------------------------------
//Twindow
//-------------------------------------------------------
TWindow::TWindow(Constants* SharedConstants){
	sharedConstants = SharedConstants;
	start = -1;
	end = -1;
	length = -1;
	sites = NULL;
	sitesInitialized = false;
	coverage = -1.0;
	fractionSitesNoData = -1.0;
	fractionsitesCoverageAtLeastTwo = -1.0;
	theta = 0.0;
	LL = 0.0;
};

TWindow::TWindow(long Start, long End, Constants* SharedConstants){
	sharedConstants = SharedConstants;
	start = Start;
	end = End;
	theta = -1.0;
	LL = 0.0;
	initSites(end - start); //end NOT in window!
};

void TWindow::initSites(long newLength){
	if(sitesInitialized)
		delete[] sites;
	length = newLength;
	sites = new TSite[length];
	sitesInitialized = true;
	coverage = -1.0;
	fractionSitesNoData = -1.0;
	fractionsitesCoverageAtLeastTwo = -1.0;
}

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

bool TWindow::addFromRead(BamTools::BamAlignment & bamAlignement){
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

	//add sites
	int internalPos = bamAlignement.Position + firstPos - start;

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
		 *                3) Function add need first pos5, then pos3
		 */

		if(bamAlignement.IsReverseStrand())
			sites[internalPos].add(bamAlignement.AlignedBases.at(pos), bamAlignement.AlignedQualities.at(pos), len - pos, pos + 1, &(sharedConstants->pmd_C_T), &(sharedConstants->pmd_G_A));
		else
			sites[internalPos].add(bamAlignement.AlignedBases.at(pos), bamAlignement.AlignedQualities.at(pos), pos + 1, len - pos, &(sharedConstants->pmd_C_T), &(sharedConstants->pmd_G_A));
	}

	//return if part of the read maps to next window
	if(lastPos == len) return false;
	else return true;
}

void TWindow::fillPGenotype(double* pGenotype, double & expTheta){
	for(int i=0; i<4; ++i){
		//homozygous genotypes
		pGenotype[sharedConstants->getGenotype(i,i)] = baseFreq[i] * (expTheta + baseFreq[i] * (1.0 - expTheta));
		//heterozygous genotypes
		for(int j=i+1; j<4; ++j){
			pGenotype[sharedConstants->getGenotype(i,j)] = 2.0 * baseFreq[i] * baseFreq[j] *  (1.0 - expTheta);
		}
	}
}

void TWindow::fillP_G(double* P_G, double* pGenotype){
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


void TWindow::calcLogLikelihood(double* pGenotype){
	LL = 0.0;
	for(int i=0; i<length; ++i){
		if(sites[i].hasData)
			LL += sites[i].calculateLogLikelihood(pGenotype);
	}
}

void TWindow::runEM(TLog* logfile){
	logfile->startIndent("Performing estimation of theta:");

	//measure runtime
	struct timeval start, end;
	gettimeofday(&start, NULL);

	//estimate initial base frequencies
	//calculate per site emission probabilities
	logfile->listFlush("Calculating emission probabilities ...");
	baseFreq.clear();
	int lengthWithData = 0;
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			sites[i].calcEmissionProbabilities();
			sites[i].addToBaseFrequencies(baseFreq);
			++lengthWithData;
		}
	}
	baseFreq.normalize();
	logfile->write(" done!");
	logfile->list("Estimating initial base frequencies ... done!");
	logfile->conclude("Pi(A) = " + toString(baseFreq[0]) + ", Pi(C) = " + toString(baseFreq[1]) + ", Pi(G) = " + toString(baseFreq[2]) + ", Pi(T) = " + toString(baseFreq[3]));

	//set initial parameters
	logfile->listFlush("Estimating initial theta ...");
	findGoodStartingTheta();
	logfile->write(" done!");
	logfile->conclude("Starting EM with theta = ", theta);
	double oldTheta, expTheta, rho, mu = lengthWithData;
	double oldLL = -9e100;

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

	//start EM loop
	logfile->listFlush("Running EM to fine tune estimations ...");
	int numThetaOnlyUpdatesDone = sharedConstants->numThetaOnlyUpdates; //do regular step first
	numThetaOnlyUpdatesDone = 0;
	int totIterations = sharedConstants->numIterations * sharedConstants->numThetaOnlyUpdates;
	for(int iter = 0; iter < totIterations; ++iter){
		//a) pre-calc expTheta
		oldTheta = theta;
		expTheta = exp(-theta);

		//b) calculate	substitution probabilities
		fillPGenotype(pGenotype, expTheta);

		//do we break EM? Check LL
		if(iter > 0 && iter % sharedConstants->numThetaOnlyUpdates == 0){
			oldLL = LL;
			calcLogLikelihood(pGenotype);
			if(LL > -9e100 && (LL - oldLL) < sharedConstants->maxEpsilon) break;
		}

		//c) Calculate all genotype probabilities for all sites
		fillP_G(P_G, pGenotype);

		//d) set rho
		/* if(iter == 0){ //estimate in first iteration
			double sumPKK = 0.0;
			for(int i=0; i<4; ++i)
				sumPKK += P_G[sharedConstants->getGenotype(i, i)];
			rho = (sumPKK) / (length - sumPKK);
		} else	*/
		rho = expTheta / (1.0 - expTheta);

		//d) Find new parameter estimates using Newton-Ralphson
		if(numThetaOnlyUpdatesDone < sharedConstants->numThetaOnlyUpdates){
			//update only theta: most difficult parameter and it is much faster to update only this one alone.
			for(int n=0; n<sharedConstants->NewtonRalphsonNumIterations; ++n){
				//i) calculate F() (Note: index is zero based!)
				F(4) = lengthWithData;
				for(int k=0; k<4; ++k){
					geno = sharedConstants->getGenotype(k, k);
					F(4) -= P_G[geno] * (rho + 1.0 ) / (rho + baseFreq[k]);
				}
				//ii) fill Jacobian (Note: index is zero based!)
				Jacobian(4,4) = 0.0;
				for(int k=0; k<4; ++k){
					tmpSum = P_G[sharedConstants->getGenotype(k, k)] / ((baseFreq[k] + rho)*(baseFreq[k] + rho));
					Jacobian(4,4) += tmpSum * (1.0 - baseFreq[k]);
				}

				//iii) now estimate new parameters
				rho = rho - F(4) / Jacobian(4,4);

				//check if we break
				if(F(4) < sharedConstants->NewtonRalphsonMaxF){
					theta = -log(rho / (1.0 + rho));
					break;
				}
			}
			++numThetaOnlyUpdatesDone;
			if(theta == oldTheta) numThetaOnlyUpdatesDone = sharedConstants->numThetaOnlyUpdates;
		} else {
			numThetaOnlyUpdatesDone = 0;
			//update all parameters in EM
			for(int n=0; n<sharedConstants->NewtonRalphsonNumIterations; ++n){
				//i) calculate F (Note: index is zero based!)
				F(4) = lengthWithData;
				F(5) = 0.0;
				for(int k=0; k<4; ++k){
					geno = sharedConstants->getGenotype(k, k);
					tmpSum = 0.0;
					for(int l=0; l<4; ++l){
						if(l != k){
							tmpSum += P_G[sharedConstants->getGenotype(k, l)];
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
					tmp[k] = P_G[sharedConstants->getGenotype(k, k)] / ((baseFreq[k] + rho)*(baseFreq[k] + rho));
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
					if(maxF < sharedConstants->NewtonRalphsonMaxF){
						theta = -log(rho / (1.0 + rho));
						break;
					}
				} else {
					++failedAttempts;
					//solve did not work -> start with higher theta!
					theta = sharedConstants->initalTheta;
					for(int i=0; i<failedAttempts; ++i)
						theta *= 10.0;
					//reset others
					mu = lengthWithData;
					LL = -9e100;
					iter = 0;
					numThetaOnlyUpdatesDone = 0;
					break;
				}
			}
		}
		//For debugging
		//std::cout << std::setprecision(9) << iter << ") theta = " << theta << "\tLL = " << LL << "\teps = " << fabs(oldLL - LL) << std::endl;
	}
	logfile->write(" done!");
	logfile->conclude("theta was estimated at ", theta);
	gettimeofday(&end, NULL);
	logfile->list("Total computation time for this window was ", end.tv_sec  - start.tv_sec, "s");
}

void TWindow::findGoodStartingTheta(){
	//assumes that initial base frequencies have been estimated and site emission probs calculated!
	double pGenotype[10];
	double initTheta = sharedConstants->initalTheta;
	double oldTheta = initTheta;
	double expTheta = exp(-initTheta);

	//calc initial LL
	fillPGenotype(pGenotype, expTheta);
	calcLogLikelihood(pGenotype);

	//run iterations
	double oldLL = LL;
	double factor = sharedConstants->initThetaSearchFactor;
	int numUpdates;
	for(int i=0; i<sharedConstants->initThetaNumSearchIterations; ++i){
		//first test increase in theta
		numUpdates = -1;
		do{
			++numUpdates;
			oldLL = LL;
			oldTheta = initTheta;
			initTheta *= factor;
			expTheta = exp(-initTheta);
			fillPGenotype(pGenotype, expTheta);
			calcLogLikelihood(pGenotype);
		} while(oldLL < LL);
		if(numUpdates == 0){
			//then test decrease in theta
			initTheta = oldTheta;
			LL = oldLL;
			//maybe smaller?
			do{
				oldLL = LL;
				oldTheta = initTheta;
				initTheta /= factor;
				expTheta = exp(-initTheta);
				fillPGenotype(pGenotype, expTheta);
				calcLogLikelihood(pGenotype);
			} while(oldLL < LL);
		}
		factor = sqrt(factor);
		initTheta = oldTheta;
		LL = oldLL;
	}
	//return previous
	theta = oldTheta;
	LL = oldLL;
}

void TWindow::calcLikelihoodSurface(std::ofstream & out, std::string & chr){
	//estimate initial base frequencies
	//calculate per site emission probabilities
	for(int i=0; i<length; ++i){
		sites[i].calcEmissionProbabilities();
		sites[i].addToBaseFrequencies(baseFreq);
	}
	baseFreq.normalize();

	//prepare storage
	double pGenotype[10];

	//calculate likelihood surface
	double minLogTheta = log10(0.000001);
	double maxLogTheta = log10(100);
	int steps = 1000;
	double stepSize = (maxLogTheta - minLogTheta) / ((double) steps - 1.0);
	double expTheta;
	double logTheta;

	for(int i=0; i<steps; ++i){
		//calc theta and expTheta
		logTheta = minLogTheta + stepSize*i;
		theta = pow(10.0, logTheta);
		expTheta = exp(-theta);

		//calculate	substitution probabilities and Likelihood
		fillPGenotype(pGenotype, expTheta);
		calcLogLikelihood(pGenotype);

		//write results
		out << chr << "\t" << start << "\t" << end << "\t" << logTheta << "\t" << LL << "\n";
	}
}


void TWindow::writeEMResults(std::ofstream & out){
	//position
	out << start << "\t" << end;
	//coverage NOTE: assumes coverage has been calculated before...
	out << "\t" << coverage << "\t" << fractionSitesNoData << "\t" << fractionsitesCoverageAtLeastTwo;
	//estimated params
	for(int i=0; i<4; ++i)
		out << "\t" << baseFreq[i];
	out << "\t" << theta << "\t" << LL << std::endl;
}

void TWindow::printPileup(std::ofstream & out, std::string & chr){
	//calc emission probs
	for(int i=0; i<length; ++i){
		sites[i].calcEmissionProbabilities();
	}
	//print pileup
	for(int i=0; i<length; ++i){
		out << chr << "\t" << start + i << "\t" << sites[i].bases.size() << "\t" << sites[i].getBases() << "\t" << sites[i].getEmissionProbs() << "\n";
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

//-------------------------------------------------------
//TGenome
//-------------------------------------------------------
TGenome::TGenome(TLog* Logfile, TParameters & params){
	logfile = Logfile;

	//read parameters
	filename = params.getParameterString("bam");
	windowSize = params.getParameterDouble("window");
	if(windowSize < 1000) throw "Window size should be at least 1Kb!";
	maxMissing = params.getParameterDoubleWithDefault("maxMissing", 0.95);

	//post mortem damage
	sharedConstants.initPmd(params, logfile);

	//EM params
	sharedConstants.numThetaOnlyUpdates = params.getParameterIntWithDefault("iterationsThetaOnly", 10);
	sharedConstants.numIterations = params.getParameterIntWithDefault("iterations", 100);

	sharedConstants.maxEpsilon = params.getParameterDoubleWithDefault("maxEps", 0.000001);
	sharedConstants.NewtonRalphsonNumIterations = params.getParameterIntWithDefault("NRiterations", 10);
	sharedConstants.NewtonRalphsonMaxF = params.getParameterDoubleWithDefault("maxF", 0.00001);

	//params regarding initial search
	sharedConstants.initalTheta = params.getParameterDoubleWithDefault("initTheta", 0.01);
	sharedConstants.initThetaSearchFactor = params.getParameterDoubleWithDefault("initThetaSearchFactor", 100);
	sharedConstants.initThetaNumSearchIterations = params.getParameterDoubleWithDefault("initThetaNumSearchIterations", 10);

	//outputname
	outname = params.getParameterStringWithDefault("out", "");
	if(outname == ""){
		//guess from filename
		outputName = filename;
		outputName = extractBeforeLast(outputName, ".");
	}

	//initialize iterators
	chrNumber = -1;
	chrLength = -1;
	curStart = -1;
	curEnd = -1;
	oldPos = -1;

	//open BAM file
	if (!bamReader.Open(filename))
		throw "Failed to open BAM file '" + filename + "'!";

	//load index file
	if(!bamReader.LocateIndex())
		throw "No index file found for BAM file '" + filename + "'!";

	//read header
	bamHeader = bamReader.GetHeader();
	chrIterator = bamHeader.Sequences.End();

	//prepare windows
	curWindow = new TWindow(&sharedConstants);
	nextWindow = new TWindow(&sharedConstants);
};

bool TGenome::iterateChromosome(){
	if(chrIterator == bamHeader.Sequences.End()){
		chrIterator = bamHeader.Sequences.Begin();
		chrNumber = 0;
	} else {
		logfile->endIndent();
		//move to next
		++chrIterator;
		++chrNumber;
	}

	//did we reach end?
	if(chrIterator == bamHeader.Sequences.End()){
		return false;
	}

	//jump reader
	bamReader.Jump(chrNumber, 0);

	//restart windows
	chrLength = stringToLong(chrIterator->Length);
	curStart = 0;
	curEnd = 0;
	oldPos = -1;
	nextWindow->move(0, windowSize);

	//write progress
	logfile->startIndent("Parsing chromosome '" + chrIterator->Name + "':");
	return true;
}

bool TGenome::iterateWindow(){
	//move to next region
	curStart = curEnd;
	curEnd += windowSize;
	if(curEnd > chrLength) curEnd = chrLength + 1;
	if(curStart > chrLength) return false;

	//swap windows and move next
	TWindow* tmp = curWindow;
	curWindow = nextWindow;
	nextWindow = tmp;
	nextWindow->move(curEnd, curEnd + windowSize);

	return true;
};

bool TGenome::readData(){
	logfile->listFlush("Reading data on '" + chrIterator->Name + "' at [" + toString(curStart) + ", " + toString(curEnd) + ") ...");

	//measure runtime
	struct timeval start, end;
	gettimeofday(&start, NULL);

	//parse through reads
	int numReads = 0;
	while(bamReader.GetNextAlignment(bamAlignement) && bamAlignement.RefID==chrNumber){
		//std::cout << "REF ID = " << bamAlignement.RefID << "\tpos = " << bamAlignement.Position << std::endl;
		//only take those reads that pass QC
		if(!bamAlignement.IsFailedQC() && !bamAlignement.IsDuplicate() && bamAlignement.Position >= curStart){
			//check if bam file is sorted
			if(bamAlignement.Position < oldPos){
				throw "BAM file must be sorted by position!";
			}
			oldPos = bamAlignement.Position;

			//check if still within current window and add to window
			if(bamAlignement.Position >= curEnd){
				//nextWindow->addFromRead(bamAlignement);
				break;
			} else {
				++numReads;
				if(curWindow->addFromRead(bamAlignement))
					//add also to next window in case reads overhangs current window -> function returns true
					nextWindow->addFromRead(bamAlignement);

			}
		}
	}

	curWindow->calcCoverage();

	gettimeofday(&end, NULL);
	logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");
	logfile->conclude("read data from " + toString(numReads) + " reads.");
	logfile->conclude("coverage is " + toString(curWindow->coverage));
	logfile->conclude(toString(curWindow->fractionsitesCoverageAtLeastTwo) + " of all sites are covered at least twice");
	logfile->conclude(toString(curWindow->fractionSitesNoData) + " of all sites have no data");

	return true;
};

void TGenome::estimateTheta(){
	//write header
	std::ofstream out;
	out.open((outputName + "_theta_estimates.txt").c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";

	//write header
	out << std::setprecision(9) << "Chr\tstart\tend\tcoverage\tmissing\ttwoOrMore\tf(A)\tf(C)\tf(G)\tf(T)\ttheta\tLL" << std::endl;

	//iterate through windows
	while(iterateChromosome()){
		while(iterateWindow()){
			//read data for current window
			readData();

			//check if we have data -> can be extended to ensure
			if(curWindow->fractionSitesNoData > maxMissing){
				logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
			} else {
				//run EM
				curWindow->runEM(logfile);

				//save results to file
				out << chrIterator->Name << "\t";
				curWindow->writeEMResults(out);
			}
		}
	}

	//close output
	out.close();
}

void TGenome::calcLikelihoodSurfaces(){
	//write header
	std::ofstream out;
	out.open((outputName + "_LLsurface.txt").c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";

	//write header
	out << "Chr\tstart\tend\tlog10Theta\tLL" << std::endl;

	//iterate through windows
	while(iterateChromosome()){
		while(iterateWindow()){
			//read data for current window
			readData();

			//check if we have data -> can be extended to ensure
			if(curWindow->fractionSitesNoData > maxMissing){
				logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
			} else {
				//calc surface
				logfile->listFlush("Calculating likelihood surface ...");
				curWindow->calcLikelihoodSurface(out, chrIterator->Name);
				logfile->write(" done!");
			}

			break;
		}
		break;
	}

	//close output
	out.close();
}

void TGenome::printPileup(){
	//write header
	std::ofstream out;
	out.open(outputName.c_str());
	if(!out) throw "Failed to open output file '" + outputName + "'!";

	//write header
	out << "Chr\tposition\tcoverage\tbases";
	for(int i=0; i<10; ++i)
		out << "\tPem(" << sharedConstants.getGenotypeString(i) << ")";
	out << std::endl;

	//iterate through windows
	while(iterateChromosome()){
		while(iterateWindow()){
			//read data for current window
			readData();

			//print pileup
			curWindow->printPileup(out, chrIterator->Name);

			break;
		}
		break;
	}

	//close output
	out.close();
}
