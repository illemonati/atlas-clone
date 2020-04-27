/*
 * TAlleleFrequencyEstimator.cpp
 *
 *  Created on: May 13, 2019
 *      Author: wegmannd
 */

#include "TAlleleFrequencyEstimator.h"

////////////////////////////////////////////////////////////////////////////////////////////////
// TAlleleHardyWeinbergFreqEstimator                                                          //
////////////////////////////////////////////////////////////////////////////////////////////////
TAlleleFreqEstimatorHardyWeinberg::TAlleleFreqEstimatorHardyWeinberg(){
	maxIter = 1000;
};

double TAlleleFreqEstimatorHardyWeinberg::estimate(const TSampleLikelihoods* storage, const uint32_t numSamplesInPop, double epsilonF, TGlfConverter & glfConverter){
	pGenotype.set(0.5);
	double weights[3];

	//run EM
	int iter = 0;
	double epsilon = epsilonF + 1.0;
	while(iter < maxIter && epsilon > epsilonF){
		double old_f = pGenotype.f;

		//calculate sums
		double sum_1 = 0.0; double sum_2 = 0.0;
		int n = 0;
		for(uint32_t i=0; i<numSamplesInPop; i++){
			if(!storage[i].isMissing){
				if(storage[i].isHaploid){
					weights[0] = glfConverter[ storage[i].glfLikelihood_0 ] * pGenotype.oneMinusf;
					weights[1] = glfConverter[ storage[i].glfLikelihood_1 ] * pGenotype.f;
					double sum = weights[0] + weights[1];

					//add to sums
					sum_1 += weights[1] / sum;
					n += 1;
				} else {
					//calculate weights
					weights[0] = glfConverter[ storage[i].glfLikelihood_0 ] * pGenotype[0];
					weights[1] = glfConverter[ storage[i].glfLikelihood_1 ] * pGenotype[1];
					weights[2] = glfConverter[ storage[i].glfLikelihood_2 ] * pGenotype[2];
					double sum = weights[0] + weights[1] + weights[2];

					//add to sums
					sum_1 += weights[1] / sum;
					sum_2 += weights[2] / sum;
					n += 2;
				}
			}
		}

		//estimate f
		pGenotype.set((sum_1 + 2.0 * sum_2) / (double) n);

		//calculate F
		epsilon = fabs(pGenotype.f - old_f);
	}

	//return estimate
	return pGenotype.f;
};

////////////////////////////////////////////////////////////////////////////////////////////////
// TAlleleHardyWeinbergFreqEstimator                                                          //
////////////////////////////////////////////////////////////////////////////////////////////////
TAlleleFreqEstimatorBayes::TAlleleFreqEstimatorBayes(TParameters & Parameters, TLog* logfile, TRandomGenerator* RandomGenerator){
	randomGenerator = RandomGenerator;
	logfile->startIndent("Initializing Bayesian allele frequency estimator:");

	//prior
	alpha = Parameters.getParameterDoubleWithDefault("alpha", 0.7);
	beta = Parameters.getParameterDoubleWithDefault("beta", 0.7);
	alphaMinusOne = alpha - 1.0;
	betaMinusOne = beta - 1.0;

	//prior support
	minPriorSupport = 1e-16;
	maxPriorSupport = 1.0 - 1e-16;
	priorDensAtMin = _prior(minPriorSupport);
	priorDensAtMax = _prior(maxPriorSupport);

	logfile->list("Will use a beta(" + toString(alpha) + "," + toString(beta) + ") prior (alpha, beta).");
	logfile->endIndent();

	//MAP estimation
	numMAPSIterations = Parameters.getParameterIntWithDefault("MAPIterations", 100);
	logfile->list("Will search for the MAP using ", numMAPSIterations, " iterations (MAPIterations).");
	f_MAP = 0.5;
	f_CI_lower = 0.0;
	f_CI_upper = 1.0;
	logDensity_atMAP = 0.0;

	//prepare initial search grid between 0.0 and 1.0
	credibleInterval = Parameters.getParameterDoubleWithDefault("credibleInterval", 0.9);
	logfile->list("Will calculate the ", credibleInterval, " Credible Interval (credibleInterval).");
	initialGridSize = Parameters.getParameterIntWithDefault("initialGridSize", 101);
	logfile->list("Will use an initial grid of size ", initialGridSize, " to identify relevant frequency range (initialGridSize).");
	if(initialGridSize < 3){
		throw "Initial grid size must be >= 3!";
	}
	initialGridLast = initialGridSize - 1;
	density_initialGrid = new double[initialGridSize];
	f_initialGrid = new double[initialGridSize];
	double step = 1.0 / (double) initialGridLast;
	f_initialGrid[0] = 0.0;
	for(int i=1; i<initialGridLast; ++i){
		f_initialGrid[i] = i * step;
	}
	f_initialGrid[initialGridLast] = 1.0;

	//final grid
	gridSize = Parameters.getParameterIntWithDefault("gridSize", 1001);
	logfile->list("Will use a grid of size ", gridSize, " to calculate credible interval (gridSize).");
	if(gridSize < 10){
		throw "Initial grid size must be >= 10!";
	}
	gridLast = gridSize - 1;
	logGridThreshold = Parameters.getParameterDoubleWithDefault("logGridThreshold", 14.0);
	logfile->list("Will use a threshold ", logGridThreshold, " to span the grid (logGridThreshold).");
	if(logGridThreshold < 1.0){
		throw "grid threshold must be >= 1.0!";
	}
	density_grid = new double[gridSize];
	f_grid = new double[gridSize];
};

TAlleleFreqEstimatorBayes::~TAlleleFreqEstimatorBayes(){
	delete[] f_initialGrid;
	delete[] density_initialGrid;
	delete[] f_grid;
	delete[] density_grid;
	delete mcmcSamples;
};

double TAlleleFreqEstimatorBayes::guessInitialAlleleFrequency(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, TGlfConverter & glfConverter){
	//calculate average posterior probs using a non-informative prior.
	double sum_1 = 0.0;
	double sum_2 = 0.0;
	int n = 0;
	for(uint32_t i=0; i<numSamplesInPopulation; i++){
		if(!storage[i].isMissing){
			if(storage[i].isHaploid){
				double sum = glfConverter[ storage[i].glfLikelihood_0 ] + glfConverter[ storage[i].glfLikelihood_1 ];

				//add to sums
				sum_1 += glfConverter[ storage[i].glfLikelihood_1 ] / sum;
				n += 1;
			} else {
				double sum = glfConverter[ storage[i].glfLikelihood_0 ] + glfConverter[ storage[i].glfLikelihood_1 ] + glfConverter[ storage[i].glfLikelihood_2 ];

				//add to sums
				sum_1 += glfConverter[ storage[i].glfLikelihood_1 ] / sum;
				sum_2 += glfConverter[ storage[i].glfLikelihood_2 ] / sum;
				n += 2;
			}
		}
	}

	return (sum_1 + 2.0 * sum_2) / (double) n;
};


double TAlleleFreqEstimatorBayes::_prior(const double & f){
	return alphaMinusOne * log(f) + betaMinusOne * log(1.0 - f);
};

double TAlleleFreqEstimatorBayes::_prior(const THardyWeinbergGenotypeProbabilities & pGenotype){
	if(pGenotype.f < minPriorSupport){
		return priorDensAtMin;
	} else if(pGenotype.f > maxPriorSupport){
		return priorDensAtMax;
	} else {
		return _prior(pGenotype.f);
	}
};

double TAlleleFreqEstimatorBayes::calcPosterior(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, THardyWeinbergGenotypeProbabilities & pGenotype, TGlfConverter & glfConverter){
	double LL = 0.0;

	for(uint32_t i=0; i<numSamplesInPopulation; i++){
		if(!storage[i].isMissing){
			if(storage[i].isHaploid){
				LL += log(glfConverter[ storage[i].glfLikelihood_0 ] * pGenotype.oneMinusf + glfConverter[ storage[i].glfLikelihood_1 ] * pGenotype.f);
			} else {
				LL += log(glfConverter[ storage[i].glfLikelihood_0 ] * pGenotype[0] + glfConverter[ storage[i].glfLikelihood_1 ] * pGenotype[1] + glfConverter[ storage[i].glfLikelihood_2 ] * pGenotype[2]);
			}
		}
	}

	//add prior
	//return LL + alphaMinusOne * log(pGenotype.f) + betaMinusOne * log(pGenotype.oneMinusf);
	//std::cout << pGenotype.f << "\t" << LL << "\t" << _prior(pGenotype) << "\t" << LL + _prior(pGenotype) << std::endl;
	return LL + _prior(pGenotype);
};

double TAlleleFreqEstimatorBayes::runMCMC(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, TGlfConverter & glfConverter, const uint32_t numIterations, const double frac){
	//prepare MCMC
	int numAccepted = 1;
	double proposalWidth = frac * (f_CI_upper - f_CI_lower);
	pGenotype.set(maxPriorSupport);
	double oldLL = calcPosterior(storage, numSamplesInPopulation, pGenotype, glfConverter);

	mcmcSamples = new double[numIterations];
	mcmcSamples[0] = maxPriorSupport;

	//run MCMC
	for(uint32_t i=1; i<numIterations; ++i){

		//propose new
		double newFreq = mcmcSamples[i-1] + randomGenerator->getRand() * proposalWidth - proposalWidth / 2.0;

		//mirror
		if(newFreq < minPriorSupport){
			newFreq = -newFreq;
		} else if(newFreq > maxPriorSupport){
			newFreq = 2.0 - newFreq;
		}

		//accept?
		pGenotype.set(newFreq);
		double newLL = calcPosterior(storage, numSamplesInPopulation, pGenotype, glfConverter);
		double h = newLL - oldLL;
		double r = log(randomGenerator->getRand());
		if(r < h){
			oldLL = newLL;
			mcmcSamples[i] = newFreq;
			++numAccepted;
		} else {
			mcmcSamples[i] = mcmcSamples[i-1];
		}
	}
	return (double) numAccepted / (double) numIterations;
};

void TAlleleFreqEstimatorBayes::fillInitialGrid(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, TGlfConverter & glfConverter){
	for(int i=0; i<initialGridSize; ++i){
		pGenotype.set(f_initialGrid[i]);
		density_initialGrid[i] = calcPosterior(storage, numSamplesInPopulation, pGenotype, glfConverter);
	}
};

void TAlleleFreqEstimatorBayes::estimateMAP(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, TGlfConverter & glfConverter){
	//use simple line-search to find MAP
	//initialize
	pGenotype.set(guessInitialAlleleFrequency(storage, numSamplesInPopulation, glfConverter));
	logDensity_atMAP = calcPosterior(storage, numSamplesInPopulation, pGenotype, glfConverter);
	double step = 0.05;

	//now do line search
	for(int i=0; i<numMAPSIterations; ++i){
		//propose new f
		pGenotype.set(pGenotype.f + step);

		//calc LL and switch if LL is lower
		double LL = calcPosterior(storage, numSamplesInPopulation, pGenotype, glfConverter);
		if(LL < logDensity_atMAP || pGenotype.f == 0.0 || pGenotype.f == 1.0){
			step = - step / 2.718282;
		}
		logDensity_atMAP = LL;
	}
	f_MAP = pGenotype.f;

	//assumes that initial grid was pre-calculated! to check and boundaries
	fillInitialGrid(storage, numSamplesInPopulation, glfConverter);

	//check if MAP is zero or one
	if(density_initialGrid[0] >= logDensity_atMAP){
		f_MAP = 0.0;
		logDensity_atMAP = density_initialGrid[0];
	} else if(density_initialGrid[initialGridLast] >= logDensity_atMAP){
		f_MAP = 1.0;
		logDensity_atMAP = density_initialGrid[initialGridLast];
	}
};

void TAlleleFreqEstimatorBayes::estimateCredibleIntervals(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, TGlfConverter & glfConverter){
	//use initial grid to define final grid
	//search first that is larger than LL_atMAP - logGridThreshold
	int first = 0;
	double relevantLL = logDensity_atMAP - logGridThreshold;
	while(density_initialGrid[first] < relevantLL && f_initialGrid[first] < f_MAP){
		++first;
	}
	if(first > 0) --first;

	//search last
	int last = initialGridLast;
	while(density_initialGrid[last] < relevantLL && f_initialGrid[last] > f_MAP){
		--last;
	}
	if(last < initialGridLast) ++last;

	//now prepare grid and calculate LL
	//make sure none is > LL_atMAP
	double step = (f_initialGrid[last] - f_initialGrid[first]) / (double) (gridLast);
	double integral = 0.0;
	for(int i=0; i<gridSize; ++i){
		f_grid[i] = f_initialGrid[first] + i * step;
		pGenotype.set(f_grid[i]);
		density_grid[i] = exp(calcPosterior(storage, numSamplesInPopulation, pGenotype, glfConverter) - logDensity_atMAP);
		integral += density_grid[i];
	}

	//adjust integral: remove half of first and last and multiply by step
	integral -= (density_grid[0] + density_grid[gridLast]) / 2.0; //first and last count half a step
	integral *= step;

	//find index left and right of MAP
	int left = gridLast;
	while(f_grid[left] >= f_MAP && left > 0){
		--left;
	}

	int right = 0;
	while(f_grid[right] <= f_MAP && right < gridLast){
		++right;
	}

	//add part from MAP to left and right grid points: average height as used when computing integral
	//then move left and right to next ones
	double halfStep = step / 2.0;
	double CI = (density_grid[left] + density_grid[right])  * halfStep;
	-- left;
	++right;

	//now find 90% CI by iteratively adding on left and right of MAP, depending on which has higher LL
	double relevantIntegral = credibleInterval * integral;
	while(CI < relevantIntegral){
		if(left < 0 || density_grid[left] < density_grid[right]){
			//add at right
			double add = (density_grid[right] + density_grid[right - 1]) * halfStep;
			if(CI + add > relevantIntegral){
				f_CI_lower = f_grid[left + 1];
				f_CI_upper = f_grid[right - 1] + step * (relevantIntegral - CI) / add;
				break;
			} else {
				CI += add;
				++right;
			}
		} else {
			//add at left
			double add = (density_grid[left] + density_grid[left + 1]) * halfStep;
			if(CI + add > relevantIntegral){
				f_CI_upper = f_grid[right - 1];
				f_CI_lower = f_grid[left + 1] + step * (relevantIntegral - CI) / add;
				break;
			} else {
				CI += add;
				--left;
			}
		}
	}
};

double TAlleleFreqEstimatorBayes::estimate(const TSampleLikelihoods* storage, const uint32_t numSamplesInPop, TGlfConverter & glfConverter){
	//get MAP estimate
	estimateMAP(storage, numSamplesInPop, glfConverter);

	//estimate credible interval
	estimateCredibleIntervals(storage, numSamplesInPop, glfConverter);

	//return MAP estimate
	return f_MAP;
};


////////////////////////////////////////////////////////////////////////////////////////////////
// TAlleleFreqEstimator                                                                       //
////////////////////////////////////////////////////////////////////////////////////////////////
TAlleleFreqEstimator::TAlleleFreqEstimator(TParameters & Parameters, TLog* Logfile){
	vcfRead = false;
	_numLoci = 0;
	logfile = Logfile;
};

std::vector<std::string> TAlleleFreqEstimator::writeHeaderAlleleFreq(bool writeGenoFreq, bool doBayesian, TAlleleFreqEstimatorBayes* BHWEstimator){
	std::vector<std::string> header = {"chr", "pos", "ref", "alt"};

	for(int p=0; p<samples.numPopulations(); p++){
		std::string name = samples.getPopulationName(p);
		header.push_back("numDiploid_" + name);
		header.push_back("numHaploid_" + name);

		if(writeGenoFreq){
			header.push_back("freqGenoRefRef_" + name);
			header.push_back("freqGenoRefAlt_" + name);
			header.push_back("freqGenoAltAlt_" + name);
			header.push_back("freqGenoRef_" + name);
			header.push_back("freqGenoAlt_" + name);
		}

		header.push_back("freqAltGF_" + name);
		header.push_back("freqAltHW_" + name);

		if(doBayesian){
			header.push_back("freqAltHW_MAP_" + name);
			header.push_back("freqAltHW_CI" + toString(BHWEstimator->credibleIntervalUsed()) + "_lower_" + name);
			header.push_back("freqAltHW_CI" + toString(BHWEstimator->credibleIntervalUsed()) + "_upper_" + name);
		}
	}
	return(header);
};

void TAlleleFreqEstimator::writeEstimatesOnePop(TOutputFileZipped & out, TGenotypeFrequencies & genoFrequencies, double alleleFrequency, TSampleLikelihoods* samples, int numSamples, TAlleleFreqEstimatorHardyWeinberg & MLHWEstimator, TAlleleFreqEstimatorBayes* BHWEstimator, double epsF, bool writeGenoFreq, bool doBayesian){
	//write num samples with data
	out << genoFrequencies.numDiploid();
	out << genoFrequencies.numHaploid();


	//write genotype frequency estimates
	if(writeGenoFreq){
		genoFrequencies.writeDiploidFrequencies(out);
		genoFrequencies.writeHaploidFrequencies(out);
	}

	//write frequency estimate based on genotype estimates
	out << alleleFrequency;

	//write HW estimates
	out << MLHWEstimator.estimate(samples, numSamples, epsF, glfConverter);

	//Bayesian estimation
	if(doBayesian){
		out << BHWEstimator->estimate(samples, numSamples, glfConverter);
		out << BHWEstimator->lowerCredibleInterval();
		out << BHWEstimator->upperCredibleInterval();
	}

};

void TAlleleFreqEstimator::estimateAlleleFreq(TParameters & Parameters, TRandomGenerator* randomGenerator){
	if(vcfRead)
		throw "VCF already read!";

	//read samples
	if(Parameters.parameterExists("samples"))
		samples.readSamples(Parameters.getParameterString("samples"), logfile);

	//create reader
	bool saveAlleleFrequencies = true;
	TPopulationLikelihoodReaderLocus reader(Parameters, logfile, saveAlleleFrequencies);
	reader.doEstimateGenotypeFrequencies();

	// open vcf file
	vcfFilename = Parameters.getParameterString("vcf");
	logfile->startIndent("Estimating allele population frequencies from VCF file '" + vcfFilename + "':");
	reader.openVCF(vcfFilename);

	//Match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		 samples.readSamplesFromVCFNames(reader.getSampleVCFNames());

	//initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);
	TPopulationLikehoodLocus storage(samples.numSamples());

	//create allele frequency estimators
	//1) Maximum likelihood HW estimator
	TAlleleFreqEstimatorHardyWeinberg MLHWEstimator;
	double epsF = Parameters.getParameterDoubleWithDefault("epsF", 0.0000001);

	//2) Maximum Likelihood genotype count estimator (use estimates from reader)
	reader.doEstimateGenotypeFrequencies();

	//3) Bayesian HW estimator (optional)
	TAlleleFreqEstimatorBayes* BHWEstimator;
	bool doBayesian = false;
	if(Parameters.parameterExists("doBayesian")){
		doBayesian = true;
		BHWEstimator = new TAlleleFreqEstimatorBayes(Parameters, logfile, randomGenerator);
	}

	bool writeGenoFreq = Parameters.parameterExists("writeGenoFreq");

	//output file
	std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
	std::string outputName = Parameters.getParameterStringWithDefault("out", tmp) + "_alleleFreq.txt.gz";
	logfile->list("Will write allele frequencies to file '" + outputName + "'.");
	TOutputFileZipped out(outputName);

	//write header
	out.writeHeader(writeHeaderAlleleFreq(writeGenoFreq, doBayesian, BHWEstimator));

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    while(reader.readDataFromVCF(storage, samples, glfConverter)){
    	//print SNP
 		reader.writePosition(out);

 		//write estimates based on genoFrequencies (if only 1 pop, use the one of reader)
 		if(samples.numPopulations() == 0){
 			writeEstimatesOnePop(out, *(reader.genotypeFrequencies()), reader.allelFrequency(), storage.samples(), samples.numSamples(), MLHWEstimator, BHWEstimator, epsF, writeGenoFreq, doBayesian);
 		} else {
 			TGenotypeFrequencies genoFrequencies;
 	 		for(int p=0; p<samples.numPopulations(); p++){
 	 			genoFrequencies.estimate(&storage[samples.startIndex(p)], samples.numSamplesInPop(p), glfConverter, epsF);
 	 			writeEstimatesOnePop(out, genoFrequencies, genoFrequencies.alleleFrequency, &storage[samples.startIndex(p)], samples.numSamplesInPop(p), MLHWEstimator, BHWEstimator, epsF, writeGenoFreq, doBayesian);
 	 		}
 		}

 		//end line
 		out << std::endl;

 		//update for next
 		++_numLoci;
     }

    //clean up
	vcfRead = true;
	out.close();
	if(doBayesian){
		delete BHWEstimator;
	}

    //report final status
	logfile->endIndent();
	reader.concludeFilters();
	if(reader.numAcceptedLoci() < 1)
		throw "No usable loci in VCF file '" + vcfFilename + "'!";
	logfile->endIndent();
};

std::vector<std::string> TAlleleFreqEstimator::writeHeaderAlleleFreqComparison(bool writeGenoFreq, TAlleleFreqEstimatorBayes* BHWEstimator){
	std::vector<std::string> header = {"chr", "pos", "ref", "alt"};

	for(int p=0; p<samples.numPopulations(); p++){
		std::string name = samples.getPopulationName(p);
		header.push_back("numDiploid_" + name);
		header.push_back("numHaploid_" + name);

		header.push_back("freqAltHW_MAP_" + name);
		header.push_back("freqAltHW_CI" + toString(BHWEstimator->credibleIntervalUsed()) + "_lower_" + name);
		header.push_back("freqAltHW_CI" + toString(BHWEstimator->credibleIntervalUsed()) + "_upper_" + name);
		header.push_back("MCMCAcceptanceRate_" + name);
	}

    //pairwise comparisons
	for(int p1=0; p1<(samples.numPopulations()-1); ++p1){
		for(int p2 = p1+1; p2 < samples.numPopulations(); ++p2){
			header.push_back("P(f_" + samples.getPopulationName(p1) + "<f_" + samples.getPopulationName(p2) + ")");
		}
	}
	return(header);
};



void TAlleleFreqEstimator::compareAlleleFreq(TParameters & Parameters, TRandomGenerator* randomGenerator){
	if(vcfRead)
		throw "VCF already read!";

	//read samples
	if(Parameters.parameterExists("samples"))
		samples.readSamples(Parameters.getParameterString("samples"), logfile);

	//create reader
	bool saveAlleleFrequencies = true;
	TPopulationLikelihoodReaderLocus reader(Parameters, logfile, saveAlleleFrequencies);
	reader.doEstimateGenotypeFrequencies();

	//open vcf file
	vcfFilename = Parameters.getParameterString("vcf");
	logfile->startIndent("Estimating allele population frequencies from VCF file '" + vcfFilename + "':");
	reader.openVCF(vcfFilename);

	//match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		 samples.readSamplesFromVCFNames(reader.getSampleVCFNames());

	if(samples.numPopulations() < 2){
		throw "Need to define at least 2 populations in samples file! Use 'task=alleleFreq' to estimate allele frequencies for a single population.";
	}

	//initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);
	TPopulationLikehoodLocus storage(samples.numSamples());

	//create bayesian allele frequency estimator
	TAlleleFreqEstimatorBayes* BHWEstimator;
	BHWEstimator = new TAlleleFreqEstimatorBayes(Parameters, logfile, randomGenerator);

	//genotype frequencies estimator
	TGenotypeFrequencies genoFrequencies;

	//variables for MCMC chains
	std::vector<double*> mcmcChains;
	int numIterations = Parameters.getParameterIntWithDefault("iterations", 1000);
	double frac = Parameters.getParameterDoubleWithDefault("proposalFrac", 3.0);
	if(numIterations < 1)
		throw "Cannot run MCMC for less than 1 iteration!";
	if(frac <= 0.0)
		throw "proposalFrac must be larger than 0!";
	logfile->list("Running MCMC for " + toString(numIterations) + " iterations (parameter 'iterations') with a propsal width of " + toString(frac) + " times the 90% confidence interval (parameter 'proposalFrac').");

	//output file
	std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
	std::string outputName = Parameters.getParameterStringWithDefault("out", tmp) + "_alleleFreqComparison.txt.gz";
	logfile->list("Will write allele frequencies to file '" + outputName + "'.");
	TOutputFileZipped out(outputName);

	//write header
	out.writeHeader(writeHeaderAlleleFreqComparison(Parameters.parameterExists("writeGenoFreq"), BHWEstimator));

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    while(reader.readDataFromVCF(storage, samples, glfConverter)){
    	//print SNP
 		reader.writePosition(out);
		logfile->list("Running MCMC's for " + toString(samples.numPopulations()) +" populations at locus " + reader.chr() + ":" + toString(reader.position()));

 		for(int p=0; p<samples.numPopulations(); p++){
			//write num samples with data
			genoFrequencies.estimate(&storage[samples.startIndex(p)], samples.numSamplesInPop(p), glfConverter, Parameters.getParameterDoubleWithDefault("epsF", 0.0000001));
			out << genoFrequencies.numDiploid();
			out << genoFrequencies.numHaploid();

			//Bayesian estimation
			out << BHWEstimator->estimate(&storage[samples.startIndex(p)], samples.numSamplesInPop(p), glfConverter);
			out << BHWEstimator->lowerCredibleInterval();
			out << BHWEstimator->upperCredibleInterval();

			//store mcmc chain
			double acceptanceRate = BHWEstimator->runMCMC(&storage[samples.startIndex(p)], samples.numSamplesInPop(p), glfConverter, numIterations, frac);
			mcmcChains.push_back(BHWEstimator->getMcmcSamples());
			out << acceptanceRate;
 		}

 		//update for next
 		++_numLoci;
     }

    //end VCF file and report final status
	vcfRead = true;
	logfile->endIndent();
	reader.concludeFilters();
	if(reader.numAcceptedLoci() < 1)
		throw "No usable loci in VCF file '" + vcfFilename + "'!";
	logfile->endIndent();

    //pairwise comparisons
	for(int p1=0; p1<(samples.numPopulations()-1); ++p1){
		//compare populations
		for(int p2 = p1+1; p2 < samples.numPopulations(); ++p2){
			logfile->list("Comparing allele frequencies of '" + samples.getPopulationName(p1) + "' and '" + samples.getPopulationName(p2)+ "'");

			int smallerThan = 0;
			for(int i=0; i<numIterations; ++i){
				if(mcmcChains[p1][i] < mcmcChains[p2][i]){
					++smallerThan;
				}
			}
			out << (double) smallerThan / (double) numIterations;
		}
	}

	//write MCMC to file?
	std::string writePop = Parameters.getParameterStringWithDefault("writeMCMC", "");
	std::vector<std::string> popToWriteMCMC;
	std::map<std::string, int> writePopMap;
	TOutputFileZipped outT;

	if(Parameters.parameterExists("writeMCMC")){
		std::string tempOut = tmp + "_MCMC.txt.gz";
		logfile->list("Will write MCMC chains to file '" + tempOut + "'.");
		TOutputFileZipped outT(tempOut);
		fillVectorFromString(writePop, popToWriteMCMC, ',');

		//write header
		outT.writeHeader(popToWriteMCMC);

		//fill map
		for(unsigned int i=0; i<popToWriteMCMC.size(); ++i){
			if(samples.populationExists(popToWriteMCMC[i]))
				writePopMap.emplace(popToWriteMCMC[i], i);
			else
				throw "Population '" + popToWriteMCMC[i] + "' does not exist!";
		}


		//write to file
		for(int i=0; i<numIterations; ++i){
			for(int p1=0; p1<(samples.numPopulations()); ++p1){
				if(writePopMap.find(samples.getPopulationName(p1)) != writePopMap.end())
					outT << mcmcChains[p1][i];
			}
			outT.endLine();
		}
		outT.close();
	}

	//clean up
	out << std::endl;
	out.close();
	delete BHWEstimator;
};

