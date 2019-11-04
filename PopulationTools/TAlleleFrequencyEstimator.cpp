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

double TAlleleFreqEstimatorHardyWeinberg::estimate(TPopulationLikehoodLocus & storage, double epsilonF, TGlfConverter & glfConverter){
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
		for(int i=0; i<storage.numSamples; i++){
			if(!storage[i].isMissing){
				if(storage[i].isHaploid){
					weights[0] = glfConverter[ storage[i][0] ] * pGenotype.oneMinusf;
					weights[1] = glfConverter[ storage[i][1] ] * pGenotype.f;
					double sum = weights[0] + weights[1];

					//add to sums
					sum_1 += weights[1] / sum;
					n += 1;
				} else {
					//calculate weights
					weights[0] = glfConverter[ storage[i][0] ] * pGenotype[0];
					weights[1] = glfConverter[ storage[i][1] ] * pGenotype[1];
					weights[2] = glfConverter[ storage[i][2] ] * pGenotype[2];
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
	alpha = Parameters.getParameterDoubleWithDefault("alpha", 0.5);
	beta = Parameters.getParameterDoubleWithDefault("beta", 0.5);
	alphaMinusOne = alpha - 1.0;
	betaMinusOne = beta - 1.0;
	logfile->list("Will use a beta(" + toString(alpha) + "," + toString(beta) + ") prior (alpha, beta).");
	logfile->endIndent();

	//MAP estimation
	numMAPSIterations = Parameters.getParameterIntWithDefault("MAPIterations", 100);
	logfile->list("Will search for the MAP using ", numMAPSIterations, " iterations (MAPIterations).");
	f_MAP = 0.5;
	f_CI_lower = 0.0;
	f_CI_upper = 1.0;
	LL_atMAP = 0.0;

	//prepare initial search grid between 0.0 and 1.0
	credibleInterval = Parameters.getParameterDoubleWithDefault("credibleInterval", 0.9);
	logfile->list("Will calculate the ", credibleInterval, " Credible Interval (credibleInterval).");
	initialGridSize = Parameters.getParameterIntWithDefault("initialGridSize", 101);
	logfile->list("Will use an initial grid of size ", initialGridSize, " to identify relevant frequency range (initialGridSize).");
	if(initialGridSize < 3){
		throw "Initial grid size must be >= 3!";
	}
	initialGridLast = initialGridSize - 1;
	LL_initialGrid = new double[initialGridSize];
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
	Likelihood_grid = new double[gridSize];
	f_grid = new double[gridSize];
};

TAlleleFreqEstimatorBayes::~TAlleleFreqEstimatorBayes(){
	delete[] f_initialGrid;
	delete[] LL_initialGrid;
	delete[] f_grid;
	delete[] Likelihood_grid;
};

double TAlleleFreqEstimatorBayes::guessInitialAlleleFrequency(TPopulationLikehoodLocus & storage, TGlfConverter & glfConverter){
	//calculate average posterior probs using a non-informative prior.
	double sum_1 = 0.0;
	double sum_2 = 0.0;
	int n = 0;
	for(int i=0; i<storage.numSamples; i++){
		if(!storage[i].isMissing){
			if(storage[i].isHaploid){
				double sum = glfConverter[ storage[i][0] ] + glfConverter[ storage[i][1] ];

				//add to sums
				sum_1 += glfConverter[ storage[i][1] ] / sum;
				n += 1;
			} else {
				double sum = glfConverter[ storage[i][0] ] + glfConverter[ storage[i][1] ] + glfConverter[ storage[i][2] ];

				//add to sums
				sum_1 += glfConverter[ storage[i][1] ] / sum;
				sum_2 += glfConverter[ storage[i][2] ] / sum;
				n += 2;
			}
		}
	}

	return (sum_1 + 2.0 * sum_2) / (double) n;
};

double TAlleleFreqEstimatorBayes::calcLL(TPopulationLikehoodLocus & storage, THardyWeinbergGenotypeProbabilities & pGenotype, TGlfConverter & glfConverter){
	double LL = 0.0;

	for(int i=0; i<storage.numSamples; i++){
		if(!storage[i].isMissing){
			if(storage[i].isHaploid){
				LL += log(glfConverter[ storage[i][0] ] * pGenotype.oneMinusf + glfConverter[ storage[i][1] ] * pGenotype.f);
			} else {
				LL += log(glfConverter[ storage[i][0] ] * pGenotype[0] + glfConverter[ storage[i][1] ] * pGenotype[1] + glfConverter[ storage[i][2] ] * pGenotype[2]);
			}
		}
	}
	return LL;
};

void TAlleleFreqEstimatorBayes::fillInitialGrid(TPopulationLikehoodLocus & storage, TGlfConverter & glfConverter){
	for(int i=0; i<initialGridSize; ++i){
		pGenotype.set(f_initialGrid[i]);
		LL_initialGrid[i] = calcLL(storage, pGenotype, glfConverter);
	}
};

void TAlleleFreqEstimatorBayes::estimateMAP(TPopulationLikehoodLocus & storage, TGlfConverter & glfConverter){
	//use simple line-search to find MAP
	//initialize
	pGenotype.set(guessInitialAlleleFrequency(storage, glfConverter));
	LL_atMAP = calcLL(storage, pGenotype, glfConverter);
	double step = 0.05;

	//now do line search
	for(int i=0; i<numMAPSIterations; ++i){
		//propose new f
		pGenotype.set(pGenotype.f + step);

		//calc LL and switch if LL is lower
		double LL = calcLL(storage, pGenotype, glfConverter);
		if(LL < LL_atMAP || pGenotype.f == 0.0 || pGenotype.f == 1.0){
			step = - step / 2.718282;
		}
		LL_atMAP = LL;
	}
	f_MAP = pGenotype.f;

	//assumes that initial grid was pre-calculated! to check and boundaries
	fillInitialGrid(storage, glfConverter);

	//check if MAP is zero or one
	if(LL_initialGrid[0] >= LL_atMAP){
		f_MAP = 0.0;
		LL_atMAP = LL_initialGrid[0];
	} else if(LL_initialGrid[initialGridLast] >= LL_atMAP){
		f_MAP = 1.0;
		LL_atMAP = LL_initialGrid[initialGridLast];
	}
};

void TAlleleFreqEstimatorBayes::estimateCredibleIntervals(TPopulationLikehoodLocus & storage, TGlfConverter & glfConverter){
	//use initial grid to define final grid
	//search first that is larger than LL_atMAP - logGridThreshold
	int first = 0;
	double relevantLL = LL_atMAP - logGridThreshold;
	while(LL_initialGrid[first] < relevantLL && f_initialGrid[first] < f_MAP){
		++first;
	}
	if(first > 0) --first;

	//search last
	int last = initialGridLast;
	while(LL_initialGrid[last] < relevantLL && f_initialGrid[last] > f_MAP){
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
		Likelihood_grid[i] = exp(calcLL(storage, pGenotype, glfConverter) - LL_atMAP);
		integral += Likelihood_grid[i];
	}

	//adjust integral: remove half of first and last and multiply by step
	integral -= (Likelihood_grid[0] + Likelihood_grid[gridLast]) / 2.0; //first and last count half a step
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
	double CI = (Likelihood_grid[left] + Likelihood_grid[right])  * halfStep;
	-- left;
	++right;

	//now find 90% CI by iteratively adding on left and right of MAP, depending on which has higher LL
	double relevantIntegral = credibleInterval * integral;
	while(CI < relevantIntegral){
		if(left < 0 || Likelihood_grid[left] < Likelihood_grid[right]){
			//add at right
			double add = (Likelihood_grid[right] + Likelihood_grid[right - 1]) * halfStep;
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
			double add = (Likelihood_grid[left] + Likelihood_grid[left + 1]) * halfStep;
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

double TAlleleFreqEstimatorBayes::estimate(TPopulationLikehoodLocus & storage, TGlfConverter & glfConverter){
	//get MAP estimate
	estimateMAP(storage, glfConverter);

	//estimate credible interval
	estimateCredibleIntervals(storage, glfConverter);

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

void TAlleleFreqEstimator::estimateAlleleFreq(TParameters & Parameters, TRandomGenerator* randomGenerator){
	if(vcfRead)
		throw "VCF already read!";

	//read samples
	TPopulationSamples samples;
	if(Parameters.parameterExists("samples"))
		samples.readSamples(Parameters.getParameterString("samples"), logfile);

	//create reader
	bool saveAlleleFrequencies = true;
	TPopulationLikelihoodReader reader(Parameters, logfile, saveAlleleFrequencies);
	reader.doEstimateGenotypeFrequencies();

	// open vcf file
	vcfFilename = Parameters.getParameterString("vcf");
	logfile->startIndent("Estimating allele population frequencies from VCF file '" + vcfFilename + "':");
	reader.openVCF(vcfFilename, logfile);

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

	//output file
	std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
	std::string outputName = Parameters.getParameterStringWithDefault("out", tmp) + "_alleleFreq.txt.gz";
	logfile->list("Will write allele frequencies to file '" + outputName + "'.");
	TOutputFileZipped out(outputName);

	//write header
	std::vector<std::string> header = {"chr", "pos", "ref", "alt", "numDiploid", "numHaploid", "freqAltHW", "freqGenoRefRef", "freqGenoRefAlt", "freqGenoAltAlt", "freqGenoRef", "freqGenoAlt", "freqAltGF"};
	if(doBayesian){
		header.push_back("freqAltHW_MAP");
		header.push_back("freqAltHW_CI" + toString(BHWEstimator->credibleIntervalUsed()) + "_lower");
		header.push_back("freqAltHW_CI" + toString(BHWEstimator->credibleIntervalUsed()) + "_upper");
	}
	out.writeHeader(header);

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    while(reader.readDataFromVCF(storage, samples, glfConverter, logfile)){
    	//print SNP
 		reader.writePosition(out);

 		//write num samples with data
 		out << reader.genotypeFrequencies()->numDiploid();
 		out << reader.genotypeFrequencies()->numHaploid();

 		//write HW estimates
 		out << MLHWEstimator.estimate(storage, epsF, glfConverter);

 		//write genotype frequency estimates
 		reader.genotypeFrequencies()->writeDiploidFrequencies(out);
 		reader.genotypeFrequencies()->writeHaploidFrequencies(out);

 		//write frequency estimate based on genotype estimates
 		out << reader.allelFrequency();

 		//Bayesian estimation
 		if(doBayesian){
 			out << BHWEstimator->estimate(storage, glfConverter);
 			out << BHWEstimator->lowerCredibleInterval();
 			out << BHWEstimator->upperCredibleInterval();
 		};

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
	reader.concludeFilters(logfile);
	if(reader.numAcceptedLoci() < 1)
		throw "No usable loci in VCF file '" + vcfFilename + "'!";
	logfile->endIndent();
};

