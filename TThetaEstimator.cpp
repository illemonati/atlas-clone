/*
 * TThetaEstimator.cpp
 *
 *  Created on: Sep 25, 2017
 *      Author: phaentu
 */

#include "TThetaEstimator.h"

//-------------------------------------------------------
//TThetaEstimatorTemporaryFile
//-------------------------------------------------------
TThetaEstimatorTemporaryFile::TThetaEstimatorTemporaryFile(std::string Filename, int numGenotypes){
	filename = Filename;
	fp = NULL;
	sizeOfData = sizeof(double) * numGenotypes;

	isOpen = false;
	isOpenForReading = false;
	isOpenForWriting = false;
	wasWritten = false;
};

void TThetaEstimatorTemporaryFile::openForWriting(){
	//if file was written, remove it
	clean();

	//now open
	fp = gzopen(filename.c_str(),"wb");
	isOpenForWriting = true;
	isOpenForReading = false;
	wasWritten = true;

	isOpen = true;
};

void TThetaEstimatorTemporaryFile::openForReading(){
	if(!wasWritten)
		throw "Can not parse temporary file: file was never written!";

	//make sure file is closed
	close();

	//now open
	fp = gzopen(filename.c_str(),"rb");
	isOpenForWriting = false;
	isOpenForReading = true;

	isOpen = true;
};

void TThetaEstimatorTemporaryFile::close(){
	if(isOpen){
		gzclose(fp);
		isOpen = false;
		isOpenForWriting = false;
		isOpenForReading = false;
	}
};

void TThetaEstimatorTemporaryFile::clean(){
	close();
	if(wasWritten){
		remove(filename.c_str());
		wasWritten = false;
	}
};

bool TThetaEstimatorTemporaryFile::isEOF(){
	if(!isOpenForReading) return true;
	return gzeof(fp);
}

void TThetaEstimatorTemporaryFile::save(double* data){
	if(!isOpenForWriting) throw "Can not add data to '" + filename + "': file is closed!";

	gzwrite(fp, data, sizeOfData);
};

bool TThetaEstimatorTemporaryFile::read(double* data){
	if(!isOpenForReading) throw "Can not read data from '" + filename + "': file is closed!";
	if(gzread(fp, data, sizeOfData)!=sizeOfData){
		//is end-of-file?
		if(gzeof(fp)) return false;

		//is error
		throw "Failed to read data from temporary file!";
	}
	return true;
};

//-------------------------------------------------------
//TThetaEstimatorData
//-------------------------------------------------------
TThetaEstimatorData::TThetaEstimatorData(int NumGenotypes){
	numGenotypes = NumGenotypes;
	pointerToData = NULL;
	numSitesCoveredTwiceOrMore = 0;
	totNumSitesAdded = 0;
	numSitesWithData = 0;
	cumulativeDepth = 0.0;

	isBootstrapped = false;
	numBootstrappedSites = false;
};

void TThetaEstimatorData::clear(){
	initialBaseFreq.clear();
	numSitesCoveredTwiceOrMore = 0;
	totNumSitesAdded = 0;
	numSitesWithData = 0;
	cumulativeDepth = 0.0;
	emptyStorage();
	clearBootstrap();
};

void TThetaEstimatorData::add(TSite & site){
	//assumes that emission probabilities were calculated!!
	++totNumSitesAdded;
	cumulativeDepth += site.depth();

	//add if site has data
	if(site.hasData){
		++numSitesWithData;
		saveSite(site);

		//add site to base frequency estimation
		site.addToBaseFrequencies(initialBaseFreq);

		//count sites covered >=2
		if(site.depth() > 1)
			++numSitesCoveredTwiceOrMore;
	}
};

void TThetaEstimatorData::fillBaseFreq(double* baseFreq){
	initialBaseFreq.normalize();
	for(int i=0; i<4; ++i)
		baseFreq[i] = initialBaseFreq[i];
};


void TThetaEstimatorData::writeHeader(std::ofstream & out){
	out << "\tdepth\tfracMissing\tfracTwoOrMore";
};

void TThetaEstimatorData::writeSize(std::ofstream & out){
	if(isBootstrapped){
		out << "\tNA\t" << (double) (totNumSitesAdded - numBootstrappedSites) / (double) totNumSitesAdded << "\tNA";	//estimated params
			out << "\t" << "NA";
	} else {
		out << "\t" << cumulativeDepth / (double) totNumSitesAdded << "\t" << (double) (totNumSitesAdded - numSitesWithData) / (double) totNumSitesAdded << "\t" << (double) numSitesCoveredTwiceOrMore / (double) totNumSitesAdded;	//estimated params
	}
};


//-------------------------------------------------------
//TThetaEstimatorDataVector
//-------------------------------------------------------
TThetaEstimatorDataVector::TThetaEstimatorDataVector(int NumGenotypes):TThetaEstimatorData(NumGenotypes){
	useTheseSites = &sites;
};

void TThetaEstimatorDataVector::saveSite(TSite & site){
	//store emission probabilities
	sites.push_back(new double[10]);
	pointerToData = *sites.rbegin();

	for(int g=0; g<numGenotypes; ++g)
		pointerToData[g] = site.emissionProbabilities[g];
};


void TThetaEstimatorDataVector::emptyStorage(){
	for(siteIt=sites.begin(); siteIt != sites.end(); ++siteIt)
		delete[] (*siteIt);
	sites.clear();
};


bool TThetaEstimatorDataVector::begin(){
	siteIt = useTheseSites->begin();
	return siteIt != useTheseSites->end();
};

bool TThetaEstimatorDataVector::next(){
	++siteIt;
	return siteIt != useTheseSites->end();
};

bool TThetaEstimatorDataVector::isEnd(){
	return siteIt == useTheseSites->end();
};

double* TThetaEstimatorDataVector::curGenotypeLikelihoods(){
	return *siteIt;
}

void TThetaEstimatorDataVector::bootstrap(TRandomGenerator & randomGenerator){
	//prepare variables
	clearBootstrap();
	int r = 0;

	//how many sites with data will we have? Draw from binomial
	double probHasData = (double) numSitesWithData / (double) totNumSitesAdded;
	numBootstrappedSites = randomGenerator.getBiomialRand(probHasData, totNumSitesAdded);

	//now pick among sites with data with replacement
	for(long i=0; i<numBootstrappedSites; ++i){
		r = randomGenerator.pickOne(numSitesWithData);
		bootstrappedSites.push_back(sites.at(r));
	}

	//set pointer
	useTheseSites = &bootstrappedSites;
	isBootstrapped = true;
};


void TThetaEstimatorDataVector::clearBootstrap(){
	useTheseSites = &sites;
	bootstrappedSites.clear();
	numBootstrappedSites = 0;
	isBootstrapped = false;
};


//-------------------------------------------------------
//TThetaEstimatorDataFile
//-------------------------------------------------------
TThetaEstimatorDataFile::TThetaEstimatorDataFile(int NumGenotypes, std::string TmpFileName):TThetaEstimatorData(NumGenotypes){
	dataFileName = TmpFileName;
	sites = new TThetaEstimatorTemporaryFile(dataFileName, numGenotypes);
	sites->openForWriting();

	pointerToData = new double[numGenotypes];

	bootstrapFileName = "bootstrap_" + dataFileName;
	bootstrappedSites = new TThetaEstimatorTemporaryFile(bootstrapFileName, numGenotypes);
	maxKforPoissonPlusOne = 100; //maximum number of bootstrapping copies to consider.
	poissonProb = new double[maxKforPoissonPlusOne];

	useTheseSites = sites;
};

void TThetaEstimatorDataFile::emptyStorage(){
	sites->clean();
};

void TThetaEstimatorDataFile::saveSite(TSite & site){
	useTheseSites->save(site.emissionProbabilities);
};

bool TThetaEstimatorDataFile::begin(){
	useTheseSites->openForReading();

	return(useTheseSites->read(pointerToData));
};

bool TThetaEstimatorDataFile::next(){
	return(useTheseSites->read(pointerToData));
};

bool TThetaEstimatorDataFile::isEnd(){
	return useTheseSites->isEOF();
};

double* TThetaEstimatorDataFile::curGenotypeLikelihoods(){
	return pointerToData;
};

void TThetaEstimatorDataFile::fillPoissonForBootstrap(const long toPick, const long remaining){
	double lambda = toPick / remaining;

	if(lambda < 1.0){
		double tmp = exp(-lambda);
		poissonProb[0] = tmp;
		for(int k=1; k<maxKforPoissonPlusOne; ++k){
			tmp = tmp * lambda / (double) k;
			poissonProb[k] = poissonProb[k-1] + tmp;
		}
	} else {
		//we need to accept!
		int n = round(lambda);
		for(int k=0; k<n; ++k)
			poissonProb[k] = 0.0;

		for(int k=n; k<maxKforPoissonPlusOne; ++k)
			poissonProb[k] = 1.0;
	}
}

void TThetaEstimatorDataFile::bootstrap(TRandomGenerator & randomGenerator){
	//make sure we start empty
	clearBootstrap();

	//how many sites with data will we have? Draw from binomial
	double probHasData = (double) numSitesWithData / (double) totNumSitesAdded;
	numBootstrappedSites = randomGenerator.getBiomialRand(probHasData, totNumSitesAdded);

	//open temporary file for bootstrap
	bootstrappedSites->openForWriting();

	//now pick among sites with data with replacement
	//For this wee need to go through the temporary file and decide for each site how many times it will be selected.
	//This is given by a Poisson distribution where lambda = n / N
	//where n = is the number of bootstrapped sites that still need to be picked and N the number of remaining sites in the file.

	//variables
	int r,i;
	long toPick = numBootstrappedSites;

	//and go...
	sites->openForReading();
	for(long l=0; l<numSitesWithData; ++l){
		sites->read(pointerToData);

		//calculate bootstrapping probabilities
		fillPoissonForBootstrap(toPick, numSitesWithData - l);

		//do we use this site in the bootstrap?
		r = randomGenerator.pickOne(maxKforPoissonPlusOne, poissonProb);
		if(r > 0){
			//save site
			for(i=0; i<r; ++i)
				bootstrappedSites->save(pointerToData);

			toPick -= r;
		}
	}

	//set pointer
	useTheseSites = bootstrappedSites;
	isBootstrapped = true;
};


void TThetaEstimatorDataFile::clearBootstrap(){
	useTheseSites = sites;
	bootstrappedSites->clean();

	numBootstrappedSites = 0;
	isBootstrapped = false;
};


//-------------------------------------------------------
//TThetaEstimator
//-------------------------------------------------------
TThetaEstimator::TThetaEstimator(TParameters & params, TLog* Logfile){
	logfile = Logfile;
	numGenotypes = 10;
	P_g_oneSite = new double[numGenotypes];

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

	//storing data in file?
	if(params.parameterExists("useTmpFile")){
		std::string dataFileName = params.getParameterStringWithDefault("useTmpFile", "temporaryDataForThetaEstimation.tmp.gz");
		logfile->list("Will write temporary data to file '" + dataFileName + "'.");
		data = new TThetaEstimatorDataFile(numGenotypes, dataFileName);
	} else data = new TThetaEstimatorDataVector(numGenotypes);

	//extra verbosity
	extraVerbose = params.parameterExists("extraVerbose");

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
	pGenotype = new double[numGenotypes];
	P_G = new double[numGenotypes];
	baseFreq = new double[4];

	//tmp stuff
	g = 0;
	sum = 0.0;
}

void TThetaEstimator::add(TSite & site){
	data->add(site);
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

void TThetaEstimator::fillP_G(){
	//assumes that pGenotype is set!
	for(g=0; g<numGenotypes; ++g)
		P_G[g] = 0.0;

	//calculate P_g for each site
	data->begin();
	do{
		sum = 0.0;
		for(g=0; g<numGenotypes; ++g){
			P_g_oneSite[g] =  data->curGenotypeLikelihoods()[g] * pGenotype[g];
			sum += P_g_oneSite[g];
		}
		for(g=0; g<numGenotypes; ++g)
			P_G[g] += P_g_oneSite[g] / sum;
	} while(data->next());
};


double TThetaEstimator::calcLogLikelihood(){
	double LL = 0.0;
	data->begin();
	do{
		sum = 0.0;
		for(g=0; g<numGenotypes; ++g)
			sum +=  data->curGenotypeLikelihoods()[g] * pGenotype[g];
		LL += log(sum);
	} while(data->next());

	return LL;
};

double TThetaEstimator::calcFisherInfo(double* deriv_pGenotype){
//sum Ri over all sites
	double FisherInfo = 0.0;
	double Ri, Ri_a, Ri_b;

	data->begin();
	do{
		//calc Ri
		Ri_a = 0.0; Ri_b = 0.0;
		for(g=0; g<numGenotypes; ++g){
			Ri_a += data->curGenotypeLikelihoods()[g] * deriv_pGenotype[g];
			Ri_b += data->curGenotypeLikelihoods()[g] * pGenotype[g];
		}
		Ri = Ri_a / Ri_b;

		//add to Fisher Info
		FisherInfo += Ri * (Ri + 1.0);
	} while(data->next());

	return(FisherInfo);
};

void TThetaEstimator::findGoodStartingTheta(){
	//set base frequencies to initial base frequencies
	data->fillBaseFreq(baseFreq);

	//variables
	double initTheta = initalTheta;
	double oldTheta = initTheta;
	double expTheta = exp(-initTheta);

	//calc initial LL
	fillPGenotype(expTheta);
	theta.LL = calcLogLikelihood();

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
	if(theta.theta < 0.1/(double) data->size()){
		theta.setTheta(0.1/(double) data->size());
	} else if(theta.theta > 1.0){
		theta.setTheta(1.0);
	}
}

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
		fillPGenotype(theta.expTheta);

		//c) Calculate all genotype probabilities for all sites
		fillP_G();

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
			theta.LL = calcLogLikelihood();
			if(theta.LL > -9e100 && (theta.LL - oldLL) < maxEpsilon) break;

			//maybe theta = 0?
			if(theta.theta < 0.1/(double) data->size()){
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

		if(extraVerbose) logfile->write("\n" + toString(iter) + ") current theta = " + toString(theta.theta));
		//For debugging
		//std::cout << std::setprecision(9) << iter << ") theta = " << thetaContainer.theta << "\tLL = " << thetaContainer.LL << "\teps = " << fabs(oldLL - thetaContainer.LL) << std::endl;
	}
	if(extraVerbose) logfile->write("EM converged, current theta = " + toString(theta.theta));

}

void TThetaEstimator::estimateConfidenceInterval(){
	//we estimate an approximate confidence interval for theta using the Fisher information
	//This function assumes that EM has already been run!

	//calculate P(g|theta, pi)
	fillPGenotype(theta.expTheta);

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

	double FisherInfo = calcFisherInfo(deriv_pGenotype);

	//estimate confidence interval
	//TODO: Fisher Info can be negative -> SQRT will be nan!
	theta.thetaConfidence = 1.96 / sqrt(FisherInfo);
}

//------------------------------------------------------------
//Functions to run estimation-
//------------------------------------------------------------
bool TThetaEstimator::estimateTheta(){
	if(data->sizeWithData() < 100){
		logfile->write("Can not estimate theta, less than 100 sites with data in this region!");
		return false;
	}

	//estimate starting parameters
	logfile->listFlush("Estimating initial parameters ...");
	findGoodStartingTheta();
	logfile->done();
	logfile->conclude("Initial base frequencies: Pi(A) = " + toString(baseFreq[0]) + ", Pi(C) = " + toString(baseFreq[1]) + ", Pi(G) = " + toString(baseFreq[2]) + ", Pi(T) = " + toString(baseFreq[3]));
	logfile->conclude("Starting EM with theta = ", theta.theta);

	//Run EM
	logfile->listFlush("Running EM to find ML estimate ...");
	runEMForTheta();
	logfile->done();
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
		fillPGenotype(expTheta);
		LL = calcLogLikelihood();

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


