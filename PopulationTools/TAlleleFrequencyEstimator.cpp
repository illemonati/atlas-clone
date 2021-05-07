/*
 * TAlleleFrequencyEstimator.cpp
 *
 *  Created on: May 13, 2019
 *      Author: wegmannd
 */

#include "TAlleleFrequencyEstimator.h"

namespace PopulationTools{

double THardyWeinbergGenotypeProbabilities::calcLogLikelihood(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, TGlfConverter & glfConverter){
	double LL = 0.0;

	for(uint32_t i=0; i<numSamplesInPopulation; i++){
		if(!storage[i].isMissing){
			if(storage[i].isHaploid){
				LL += log(glfConverter[ storage[i].glfLikelihood_0 ] * oneMinusf + glfConverter[ storage[i].glfLikelihood_1 ] * f);
			} else {
				LL += log(glfConverter[ storage[i].glfLikelihood_0 ] * genotypeProbabilities[0] + glfConverter[ storage[i].glfLikelihood_1 ] * genotypeProbabilities[1] + glfConverter[ storage[i].glfLikelihood_2 ] * genotypeProbabilities[2]);
			}
		}
	}

	return(LL);
};

////////////////////////////////////////////////////////////////////////////////////////////////
// TAlleleFreqEstimatorHardyWeinberg                                                          //
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
TAlleleFreqEstimatorBayes::TAlleleFreqEstimatorBayes(TParameters & Parameters, TLog* logfile, TGlfConverter* GlfConverter, TRandomGenerator* RandomGenerator){
	logfile->startIndent("Initializing Bayesian allele frequency estimator:");

	randomGenerator = RandomGenerator;
	glfConverter = GlfConverter;

	//prior
	alpha = Parameters.getParameterWithDefault("alpha", 0.7);
	beta = Parameters.getParameterWithDefault("beta", 0.7);
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
	numMAPSIterations = Parameters.getParameterWithDefault<int>("MAPIterations", 100);
	logfile->list("Will search for the MAP using ", numMAPSIterations, " iterations (MAPIterations).");
	f_MAP = 0.5;
	f_CI_lower = 0.0;
	f_CI_upper = 1.0;
	logDensity_atMAP = 0.0;

	//prepare initial search grid between 0.0 and 1.0
	credibleInterval = Parameters.getParameterWithDefault("credibleInterval", 0.9);
	logfile->list("Will calculate the ", credibleInterval, " Credible Interval (credibleInterval).");
	initialGridSize = Parameters.getParameterWithDefault<int>("initialGridSize", 101);
	logfile->list("Will use an initial grid of size ", initialGridSize, " to identify relevant frequency range (initialGridSize).");
	if(initialGridSize < 3){
		throw "Initial grid size must be >= 3!";
	}
	initialGridLast = initialGridSize - 1;
	density_initialGrid = new double[initialGridSize];
	f_initialGrid = new double[initialGridSize];
	double step = 1.0 / (double) initialGridLast;
	f_initialGrid[0] = minPriorSupport;
	for(int i=1; i<initialGridLast; ++i){
		f_initialGrid[i] = i * step;
	}
	f_initialGrid[initialGridLast] = maxPriorSupport;

	//final grid
	gridSize = Parameters.getParameterWithDefault<int>("gridSize", 1001);
	logfile->list("Will use a grid of size ", gridSize, " to calculate credible interval (gridSize).");
	if(gridSize < 10){
		throw "Initial grid size must be >= 10!";
	}
	gridLast = gridSize - 1;
	logGridThreshold = Parameters.getParameterWithDefault("logGridThreshold", 14.0);
	logfile->list("Will use a threshold ", logGridThreshold, " to span the grid (logGridThreshold).");
	if(logGridThreshold < 1.0){
		throw "grid threshold must be >= 1.0!";
	}
	density_grid = new double[gridSize];
	f_grid = new double[gridSize];

	//mcmcSamples
	mcmcSamplesInitialized = false;
};

TAlleleFreqEstimatorBayes::~TAlleleFreqEstimatorBayes(){
	delete[] f_initialGrid;
	delete[] density_initialGrid;
	delete[] f_grid;
	delete[] density_grid;
};

double TAlleleFreqEstimatorBayes::_guessInitialAlleleFrequency(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation){
	//calculate average posterior probs using a non-informative prior.
	double sum_1 = 0.0;
	double sum_2 = 0.0;
	int n = 0;

	for(uint32_t i=0; i<numSamplesInPopulation; i++){
		if(!storage[i].isMissing){
			if(storage[i].isHaploid){
				double sum = glfConverter->toScaledLikelihood( storage[i].glfLikelihood_0 ) + glfConverter->toScaledLikelihood( storage[i].glfLikelihood_1 );

				//add to sums
				sum_1 += glfConverter->toScaledLikelihood( storage[i].glfLikelihood_1 ) / sum;
				n += 1;
			} else {
				double sum = glfConverter->toScaledLikelihood( storage[i].glfLikelihood_0 ) + glfConverter->toScaledLikelihood( storage[i].glfLikelihood_1 ) + glfConverter->toScaledLikelihood( storage[i].glfLikelihood_2 );

				//add to sums
				sum_1 += glfConverter->toScaledLikelihood( storage[i].glfLikelihood_1 ) / sum;
				sum_2 += glfConverter->toScaledLikelihood( storage[i].glfLikelihood_2 ) / sum;
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

double TAlleleFreqEstimatorBayes::_calcPosterior(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, THardyWeinbergGenotypeProbabilities & pGenotype){
	return pGenotype.calcLogLikelihood(storage, numSamplesInPopulation, *glfConverter) + _prior(pGenotype);
};

void TAlleleFreqEstimatorBayes::_fillInitialGrid(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation){
	for(int i=0; i<initialGridSize; ++i){
		pGenotype.set(f_initialGrid[i]);
		density_initialGrid[i] = _calcPosterior(storage, numSamplesInPopulation, pGenotype);
	}
};

void TAlleleFreqEstimatorBayes::_estimateMAP(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation){
	//use simple line-search to find MAP
	//initialize
	pGenotype.set(_guessInitialAlleleFrequency(storage, numSamplesInPopulation));
	logDensity_atMAP = _calcPosterior(storage, numSamplesInPopulation, pGenotype);
	double step = 0.05;

	//now do line search
	for(int i=0; i<numMAPSIterations; ++i){
		//propose new f
		pGenotype.set(pGenotype.f + step);

		//calc LL and switch if LL is lower
		double LL = _calcPosterior(storage, numSamplesInPopulation, pGenotype);
		if(LL < logDensity_atMAP || pGenotype.f == 0.0 || pGenotype.f == 1.0){
			step = - step / 2.718282;
		}
		logDensity_atMAP = LL;
	}
	f_MAP = pGenotype.f;

	//assumes that initial grid was pre-calculated! to check and boundaries
	_fillInitialGrid(storage, numSamplesInPopulation);

	//check if MAP is zero or one
	if(density_initialGrid[0] >= logDensity_atMAP){
		f_MAP = 0.0;
		logDensity_atMAP = density_initialGrid[0];
	} else if(density_initialGrid[initialGridLast] >= logDensity_atMAP){
		f_MAP = 1.0;
		logDensity_atMAP = density_initialGrid[initialGridLast];
	}
};

void TAlleleFreqEstimatorBayes::_estimateCredibleIntervals(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation){
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
		density_grid[i] = exp(_calcPosterior(storage, numSamplesInPopulation, pGenotype) - logDensity_atMAP);
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
	double CI = 0;
	if(left > 0){
		CI += density_grid[left] * halfStep;
		-- left;
	}
	if(right < gridLast){
		CI += density_grid[right] * halfStep;
		++right;
	}

	//now find 90% CI by iteratively adding on left and right of MAP, depending on which has higher LL
	double relevantIntegral = credibleInterval * integral;
	if(CI >= relevantIntegral){
		f_CI_upper = f_grid[right];
		f_CI_lower = f_grid[left];
	} else {
		while(CI < relevantIntegral){
			if(right < gridSize && (left < 0 || density_grid[left] < density_grid[right] )){
				//add at right
				double add = (density_grid[right] + density_grid[right - 1]) * halfStep;
				if(CI + add >= relevantIntegral){
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
				if(CI + add >= relevantIntegral){
					f_CI_upper = f_grid[right - 1];
					f_CI_lower = f_grid[left + 1] + step * (relevantIntegral - CI) / add;
					break;
				} else {
					CI += add;
					--left;
				}
			}
		}
	}

	if(f_CI_lower == minPriorSupport)
		f_CI_lower = 0.0;
	if(f_CI_upper == maxPriorSupport)
		f_CI_upper = 1.0;
};

double TAlleleFreqEstimatorBayes::estimate(const TSampleLikelihoods* storage, const uint32_t numSamplesInPop){
	//get MAP estimate
	_estimateMAP(storage, numSamplesInPop);

	//estimate credible interval
	_estimateCredibleIntervals(storage, numSamplesInPop);

	//return MAP estimate
	return f_MAP;
};

void TAlleleFreqEstimatorBayes::composeHeader(std::vector<std::string> & header, const std::string & popName){
	header.push_back("freqAltHW_MAP_" + popName);
	header.push_back("freqAltHW_CI" + toString(credibleInterval) + "_lower_" + popName);
	header.push_back("freqAltHW_CI" + toString(credibleInterval) + "_upper_" + popName);
};

void TAlleleFreqEstimatorBayes::estimateAndWrite(const TSampleLikelihoods* storage, const uint32_t numSamplesInPop, TOutputFile & out){
	out << estimate(storage, numSamplesInPop); //MAP
	out << lowerCredibleInterval();
	out << upperCredibleInterval();
};

double TAlleleFreqEstimatorBayes::runMCMC(const TSampleLikelihoods* storage, const uint32_t numSamplesInPopulation, const double frac, std::vector<double> & mcmcSamples){
	//prepare MCMC
	size_t numAccepted = 1;
	double proposalWidth = frac * (f_CI_upper - f_CI_lower);
	if(f_MAP == 0.0)
		f_MAP = minPriorSupport;
	if(f_MAP == 1.0)
		f_MAP = maxPriorSupport;
	pGenotype.set(f_MAP);
	double oldLL = _calcPosterior(storage, numSamplesInPopulation, pGenotype);

	//run MCMC
	mcmcSamples[0] = f_MAP;
	for(size_t i=1; i<mcmcSamples.size(); ++i){

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
		double newLL = _calcPosterior(storage, numSamplesInPopulation, pGenotype);
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
	return (double) numAccepted / (double) mcmcSamples.size();
};


double TAlleleFreqEstimatorBayes::calcPosteriorf1smallerf2(std::vector<double> & mcmc1, std::vector<double> & mcmc2){
	size_t smallerThan = 0;
	size_t len = std::min(mcmc1.size(), mcmc2.size());
	for(size_t i=0; i<len; ++i){
		if(mcmc1[i] < mcmc2[i]){
			++smallerThan;
		}
	}
	return (double) smallerThan / (double) len;
};

////////////////////////////////////////////////////////////////////////////////////////////////
// TAlleleFreqMCMCOutput                                                                       //
////////////////////////////////////////////////////////////////////////////////////////////////
void TAlleleFreqMCMCOutput::initialize(std::string popString, TPopulationSamples & samples, std::string OutputName, TLog* logfile){
	//parse string to identify pops for which MCMC shoudl be written
	std::vector<std::string> tmp;
	fillContainerFromString(popString, tmp, ',');
	outputName = OutputName;

	//extract indexes
	if(tmp.size() > 0){
		//write all?
		if(tmp.size() == 1 && tmp[0] == "all"){
			for(int p=0; p<samples.numPopulations(); p++){
				popIndex.push_back(p);
				header.push_back(samples.getPopulationName(p));
			}
		} else {
			for(auto& name : tmp){
				if(samples.populationExists(name)){
					popIndex.push_back(samples.getPopulationIndex(name));
					header.push_back(name);
				} else {
					throw "Can not write MCMC: population '" + name + "' does not exist!";
				}
			}
		}

		if(popIndex.size() > 0){
			logfile->startIndent("Will write the MCMC of the following populations to files '" + outputName + "[chr]_[locus].txt.gz':");
			for(auto& name : header){
				logfile->list(name);
			}
			logfile->endIndent();
		}
	}
};

void TAlleleFreqMCMCOutput::write(std::vector< std::vector<double> > & mcmc, const std::string chr, const uint64_t pos){
	if(popIndex.size() > 0){
		//open output file
		std::string filename = outputName + chr + "_" + toString(pos) + ".txt.gz";
		outFile.open(filename);
		outFile.writeHeader(header);

		//get min length of chains
		size_t len = mcmc[0].size();
		for(auto& it : mcmc){
			if(it.size() < len)
				len = it.size();
		}

		//write MCMC to file
		for(size_t i=0; i<len; ++i){
			for(auto& p : popIndex){
				outFile << mcmc[p][i];
			}
			outFile.endLine();
		}

		//close file
		outFile.close();
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////
// TAlleleFreqEstimator                                                                       //
////////////////////////////////////////////////////////////////////////////////////////////////
TAlleleFreqEstimator::TAlleleFreqEstimator(TParameters & Parameters, TLog* Logfile){
	vcfRead = false;
	logfile = Logfile;
};

std::vector<std::string> TAlleleFreqEstimator::_composeHeaderAlleleFreq(bool writeGenoFreq, bool doBayesian, TAlleleFreqEstimatorBayes* BHWEstimator){
	std::vector<std::string> header = {"chr", "pos", "ref", "alt"};

	for(int p=0; p<samples.numPopulations(); p++){
		std::string popName = samples.getPopulationName(p);
		header.push_back("numDiploid_" + popName);
		header.push_back("numHaploid_" + popName);

		if(writeGenoFreq){
			header.push_back("freqGenoRefRef_" + popName);
			header.push_back("freqGenoRefAlt_" + popName);
			header.push_back("freqGenoAltAlt_" + popName);
			header.push_back("freqGenoRef_" + popName);
			header.push_back("freqGenoAlt_" + popName);
		}

		header.push_back("freqAltGF_" + popName);
		header.push_back("freqAltHW_" + popName);

		if(doBayesian){
			BHWEstimator->composeHeader(header, popName);
		}
	}
	return(header);
};


void TAlleleFreqEstimator::_writeBayesianEstimatesOnePop(TOutputFile & out, TSampleLikelihoods* samples, const uint32_t numSamples, TAlleleFreqEstimatorBayes* BHWEstimator){
	out << BHWEstimator->estimate(samples, numSamples); //MAP
	out << BHWEstimator->lowerCredibleInterval();
	out << BHWEstimator->upperCredibleInterval();
};

void TAlleleFreqEstimator::_writeEstimatesOnePop(TOutputFile & out, TGenotypeFrequencies & genoFrequencies, double alleleFrequency, TSampleLikelihoods* samples, uint32_t numSamples, TAlleleFreqEstimatorHardyWeinberg & MLHWEstimator, TAlleleFreqEstimatorBayes* BHWEstimator, double epsF, bool writeGenoFreq, bool doBayesian){
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
		_writeBayesianEstimatesOnePop(out, samples, numSamples, BHWEstimator);
	}

};

void TAlleleFreqEstimator::_openVCF(TParameters & Parameters){
	if(vcfRead)
		throw "VCF already read!";

	//read samples
	if(Parameters.parameterExists("samples"))
		samples.readSamples(Parameters.getParameter<std::string>("samples"), logfile);

	//create reader
	bool saveAlleleFrequencies = true;
	reader.initialize(Parameters, logfile, saveAlleleFrequencies);
	reader.doEstimateGenotypeFrequencies();

	// open vcf file
	vcfFilename = Parameters.getParameter<std::string>("vcf");
	logfile->startIndent("Estimating allele population frequencies from VCF file '" + vcfFilename + "':");
	reader.openVCF(vcfFilename);

	//Match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		samples.readSamplesFromVCFNames(reader.getSampleVCFNames());

	//initialize variables for vcf-file
	storage.resize(samples.numSamples());
};

void TAlleleFreqEstimator::_closeVCF(){
	vcfRead = true;

	//report final status
	logfile->endIndent();
	reader.concludeFilters();
	if(reader.numAcceptedLoci() < 1)
		throw "No usable loci in VCF file '" + vcfFilename + "'!";
	logfile->endIndent();
};

void TAlleleFreqEstimator::estimateAlleleFreq(TParameters & Parameters, TRandomGenerator* randomGenerator){
	//open VCF for reading
	_openVCF(Parameters);

	//create allele frequency estimators
	//1) Maximum likelihood HW estimator
	TAlleleFreqEstimatorHardyWeinberg MLHWEstimator;
	double epsF = Parameters.getParameterWithDefault("epsF", 0.0000001);

	//2) Maximum Likelihood genotype count estimator (use estimates from reader)
	reader.doEstimateGenotypeFrequencies();

	//3) Bayesian HW estimator (optional)
	TAlleleFreqEstimatorBayes* BHWEstimator;
	bool doBayesian = false;
	if(Parameters.parameterExists("doBayesian")){
		doBayesian = true;
		BHWEstimator = new TAlleleFreqEstimatorBayes(Parameters, logfile, &glfConverter, randomGenerator);
	}

	bool writeGenoFreq = Parameters.parameterExists("writeGenoFreq");

	//output file
	std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
	std::string outputName = Parameters.getParameterWithDefault<std::string>("out", tmp) + "_alleleFreq.txt.gz";
	logfile->list("Will write allele frequencies to file '" + outputName + "'.");
	TOutputFile out(outputName);

	//write header
	out.writeHeader(_composeHeaderAlleleFreq(writeGenoFreq, doBayesian, BHWEstimator));

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    while(reader.readDataFromVCF(storage, samples, glfConverter)){
    	//print SNP
 		reader.writePosition(out);

 		//write estimates based on genoFrequencies (if only 1 pop, use the one of reader)
 		if(samples.numPopulations() == 1){
 			_writeEstimatesOnePop(out, *(reader.genotypeFrequencies()), reader.allelFrequency(), storage.samples(), samples.numSamples(), MLHWEstimator, BHWEstimator, epsF, writeGenoFreq, doBayesian);
 		} else {
 			TGenotypeFrequencies genoFrequencies;
 	 		for(int p=0; p<samples.numPopulations(); p++){
 	 			genoFrequencies.estimate(&storage[samples.startIndex(p)], samples.numSamplesInPop(p), glfConverter, epsF);
 	 			_writeEstimatesOnePop(out, genoFrequencies, genoFrequencies.alleleFrequency, &storage[samples.startIndex(p)], samples.numSamplesInPop(p), MLHWEstimator, BHWEstimator, epsF, writeGenoFreq, doBayesian);
 	 		}
 		}

 		//end line
 		out << std::endl;
     }

    //clean up
	out.close();
	if(doBayesian){
		delete BHWEstimator;
	}

	//close VCF
	_closeVCF();
};

std::vector<std::string> TAlleleFreqEstimator::_composeHeaderAlleleFreqComparison(TAlleleFreqEstimatorBayes & BHWEstimator){
	std::vector<std::string> header = {"chr", "pos", "ref", "alt"};

	for(int p=0; p<samples.numPopulations(); p++){
		std::string popName = samples.getPopulationName(p);
		header.push_back("numDiploid_" + popName);
		header.push_back("numHaploid_" + popName);

		BHWEstimator.composeHeader(header, popName);

		header.push_back("MCMCAcceptanceRate_" + popName);
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
	//open VCF for reading
	_openVCF(Parameters);
	if(samples.numPopulations() < 2){
		throw "Need to define at least 2 populations in samples file! Use 'task=alleleFreq' to estimate allele frequencies for a single population.";
	}

	//create Bayesian allele frequency estimator
	TAlleleFreqEstimatorBayes BHWEstimator(Parameters, logfile, &glfConverter, randomGenerator);

	//genotype frequencies estimator
	TGenotypeFrequencies genoFrequencies;

	//variables for MCMC chains
	int numIterations = Parameters.getParameterWithDefault<int>("iterations", 100000);
	double frac = Parameters.getParameterWithDefault("proposalFrac", 3.0);
	if(numIterations < 1)
		throw "Cannot run MCMC for less than 1 iteration!";
	if(frac <= 0.0)
		throw "proposalFrac must be larger than 0!";
	logfile->list("Running MCMC for " + toString(numIterations) + " iterations with a propsal width of " + toString(frac) + " times the 90% confidence interval.");

	//prepare MCMC storage
	std::vector< std::vector<double> > mcmcChains(samples.numPopulations());
	for(int p=0; p<samples.numPopulations(); p++){
		mcmcChains[p].resize(numIterations);
	}

	//output file
	std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
	std::string outputName = Parameters.getParameterWithDefault<std::string>("out", tmp);
	logfile->list("Will write allele frequencies to file '" + outputName  + "_alleleFreqComparison.txt.gz" + "'.");
	TOutputFile out(outputName + "_alleleFreqComparison.txt.gz");
	out.writeHeader(_composeHeaderAlleleFreqComparison(BHWEstimator));

	//write MCMC to file?
	TAlleleFreqMCMCOutput traces;
	if(Parameters.parameterExists("writeMCMC")){
		traces.initialize(Parameters.getParameter<std::string>("writeMCMC"), samples, outputName + "_alleleFreq_MCMC_", logfile);
	}

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    while(reader.readDataFromVCF(storage, samples, glfConverter)){
    	//print SNP
 		reader.writePosition(out);

 		//run MCMC
 		logfile->listFlush("Running estimates for " + reader.chr() + ":" + toString(reader.position()) + " ...");
 		for(int p=0; p<samples.numPopulations(); p++){
			//write num samples with data
			genoFrequencies.estimate(&storage[samples.startIndex(p)], samples.numSamplesInPop(p), glfConverter, Parameters.getParameterWithDefault("epsF", 0.0000001));
			out << genoFrequencies.numDiploid();
			out << genoFrequencies.numHaploid();

			//Bayesian estimation
			BHWEstimator.estimateAndWrite(&storage[samples.startIndex(p)], samples.numSamplesInPop(p), out);

			//run mcmc chain for each population
			double acceptanceRate = BHWEstimator.runMCMC(&storage[samples.startIndex(p)], samples.numSamplesInPop(p), frac, mcmcChains[p]);
			out << acceptanceRate;
			traces.write(mcmcChains, reader.chr(), reader.position());
 		}

 		//do pairwise comparisons
 		for(int p1=0; p1<(samples.numPopulations()-1); ++p1){
			for(int p2 = p1+1; p2 < samples.numPopulations(); ++p2){
				out << BHWEstimator.calcPosteriorf1smallerf2(mcmcChains[p1], mcmcChains[p2]);
			}
		}
 		out << std::endl;
 		logfile->done();
    }

	//close VCF and output file
    out.close();
    _closeVCF();
};

void TAlleleFreqEstimator::writeAlleleFrequencyLikelihoods(TParameters & Parameters, TRandomGenerator* randomGenerator){
	//calculating P(D|f) at predefined f
	//open VCF for reading
	_openVCF(Parameters);

	//get vector of allele frequencies at which to calculate likelihood
	int numFreq = Parameters.getParameterWithDefault<int>("numFreq", 101);
	logfile->list("Will calculate allele frequency likelihoods at " + toString(numFreq) + " uniformly spaced frequencies.");
	double step = 1.0 / (double) (numFreq - 1);
	std::vector<double> freq(numFreq);
	std::vector<std::string> header = {"chr", "pos"};
	for(int i=0; i<numFreq; ++i){
		freq[i] = i * step;
		header.push_back("LL_" + toString(freq[i]));
	}

	//output file
	std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
	std::string outputName = Parameters.getParameterWithDefault<std::string>("out", tmp) + "_alleleFreqLikelihoods";
	logfile->list("Will write allele frequencies to files '" + outputName + "[POP].txt.gz'.");

	std::vector<TOutputFile> out(samples.numPopulations());
	if(samples.numPopulations() == 1){
		out[0].open(outputName + ".txt.gz");
		out[0].writeHeader(header);
		logfile->list("Will write allele frequency likelihoods to file '" + out[0].name() + "'.");
	} else {
		logfile->startIndent("Will write allele frequency likelihoods to files:");
		for(int p=0; p<samples.numPopulations(); p++){
			out[p].open(outputName + "_" + samples.getPopulationName(p) + ".txt.gz");
			out[p].writeHeader(header);
			logfile->list(out[p].name());
		}
		logfile->endIndent();
	}

	//prepare genotype probability object
	THardyWeinbergGenotypeProbabilities genoProb;

	//run through VCF file
	logfile->startIndent("Parsing VCF file:");
	while(reader.readDataFromVCF(storage, samples, glfConverter)){
		//calculate and write allele frequency likelihoods for every population
		for(int p=0; p<samples.numPopulations(); p++){
			out[p] << reader.chr() << reader.position();
			for(auto& f : freq){
				genoProb.set(f);
				out[p] << genoProb.calcLogLikelihood(&storage[samples.startIndex(p)], samples.numSamplesInPop(p), glfConverter);
			}
			out[p] << std::endl;
		}
	}

	//clean up
	for(int p=0; p<samples.numPopulations(); p++){
		out[p].close();
	}
	_closeVCF();

};

}; //end namespace
