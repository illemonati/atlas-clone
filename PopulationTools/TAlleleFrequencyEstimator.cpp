/*
 * TAlleleFrequencyEstimator.cpp
 *
 *  Created on: May 13, 2019
 *      Author: wegmannd
 */

#include "TAlleleFrequencyEstimator.h"

#include <math.h>
#include <stddef.h>
#include <algorithm>
#include <cstdint>
#include <exception>
#include <memory>
#include <ostream>
#include <string>

#include "GenotypeTypes.h"
#include "TGenotypeFrequencies.h"
#include "stringFunctions.h"
#include "weakTypes.h"

namespace PopulationTools{

using coretools::TParameters;
using coretools::TLog;
using coretools::TRandomGenerator;
using coretools::Probability;
using coretools::TOutputFile;
using coretools::str::toString;

////////////////////////////////////////////////////////////////////////////////////////////////
// TAlleleFreqEstimatorHardyWeinberg                                                          //
////////////////////////////////////////////////////////////////////////////////////////////////
TAlleleFreqEstimatorHardyWeinberg::TAlleleFreqEstimatorHardyWeinberg(){
	maxIter = 1000;
};

Probability TAlleleFreqEstimatorHardyWeinberg::estimate(const TSampleLikelihoods* storage, uint32_t numSamplesInPop, double epsilonF){
	using BG = genometools::BiallelicGenotype;
	genometools::THardyWeinbergGenotypeProbabilities pGenotype;
	Probability weights[3];

	//run EM
	size_t iter = 0;
	double epsilon = epsilonF + 1.0;
	while(iter < maxIter && epsilon > epsilonF){
		Probability old_f = pGenotype.f();

		//calculate sums
		double sum_1 = 0.0; double sum_2 = 0.0;
		int n = 0;
		for(uint32_t i=0; i<numSamplesInPop; i++){
			if(!storage[i].isMissing()){
				if(storage[i].isHaploid()){
					weights[0] = (Probability) storage[i][BG::haploidFirst]  * pGenotype[BG::haploidFirst];
					weights[1] = (Probability) storage[i][BG::haploidSecond]  * pGenotype[BG::haploidSecond];
					double sum = weights[0] + weights[1];

					//add to sums
					sum_1 += weights[1] / sum;
					n += 1;
				} else {
					//calculate weights
					weights[0] = (Probability) storage[i][BG::homoFirst]  * pGenotype[BG::homoFirst];
					weights[1] = (Probability) storage[i][BG::het]  * pGenotype[BG::het];
					weights[2] = (Probability) storage[i][BG::homoSecond]  * pGenotype[BG::homoSecond];
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
		epsilon = fabs(pGenotype.f() - old_f);
	}

	//return estimate
	return pGenotype.f();
};

////////////////////////////////////////////////////////////////////////////////////////////////
// TAlleleHardyWeinbergFreqEstimator                                                          //
////////////////////////////////////////////////////////////////////////////////////////////////
TAlleleFreqEstimatorBayes::TAlleleFreqEstimatorBayes(TParameters & Parameters, TLog* logfile, TRandomGenerator* RandomGenerator){
	logfile->startIndent("Initializing Bayesian allele frequency estimator:");

	randomGenerator = RandomGenerator;

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

	logfile->list("Will use a beta(", alpha, ",", beta, ") prior (alpha, beta).");
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
	_initialGrid.resize(initialGridSize);
	double step = 1.0 / (double) initialGridLast;
	_initialGrid[0].f = minPriorSupport;
	for(int i=1; i<initialGridLast; ++i){
		_initialGrid[i].f = i * step;
	}
	_initialGrid[initialGridLast].f = maxPriorSupport;

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
	_grid.resize(gridSize);
};

Probability TAlleleFreqEstimatorBayes::_guessInitialAlleleFrequency(const TSampleLikelihoods* storage, uint32_t numSamplesInPopulation){
	using BG = genometools::BiallelicGenotype;
	//calculate average posterior probs using a non-informative prior.
	double sum_1 = 0.0;
	double sum_2 = 0.0;
	int n = 0;

	for(uint32_t i=0; i<numSamplesInPopulation; i++){
		if(!storage[i].isMissing()){
			if(storage[i].isHaploid()){
				double sum = (Probability) storage[i][BG::homoFirst] + (Probability) storage[i][BG::homoSecond];

				//add to sums
				sum_1 += (Probability) storage[i][BG::homoSecond] / sum;
				n += 1;
			} else {
				double sum = (Probability) storage[i][BG::homoFirst] + (Probability) storage[i][BG::het] + (Probability) storage[i][BG::homoSecond];

				//add to sums
				sum_1 += (Probability) storage[i][BG::het] / sum;
				sum_2 += (Probability) storage[i][BG::homoSecond] / sum;
				n += 2;
			}
		}
	}

	return (sum_1 + 2.0 * sum_2) / (double) n;
};


double TAlleleFreqEstimatorBayes::_prior(const Probability & f) const {
	return alphaMinusOne * log(f) + betaMinusOne * log(f.complement());
};

double TAlleleFreqEstimatorBayes::_prior(const genometools::THardyWeinbergGenotypeProbabilities & pGenotype) const {
	if(pGenotype.f() < minPriorSupport){
		return priorDensAtMin;
	} else if(pGenotype.f() > maxPriorSupport){
		return priorDensAtMax;
	} else {
		return _prior(pGenotype.f());
	}
};

coretools::LogProbability TAlleleFreqEstimatorBayes::_calcLogLikelihood(const TSampleLikelihoods* storage, uint32_t numSamplesInPopulation, const genometools::THardyWeinbergGenotypeProbabilities & pGenotype){
	coretools::LogProbability LL = 0.0;
	for(uint32_t i=0; i<numSamplesInPopulation; i++){
		LL += storage[i].HWESum<coretools::LogProbability>(pGenotype);
	}
	return(LL);
};

double TAlleleFreqEstimatorBayes::_calcPosterior(const TSampleLikelihoods* storage, uint32_t numSamplesInPopulation, const genometools::THardyWeinbergGenotypeProbabilities & pGenotype){
	return _calcLogLikelihood(storage, numSamplesInPopulation, pGenotype) + _prior(pGenotype);
};

void TAlleleFreqEstimatorBayes::_fillInitialGrid(const TSampleLikelihoods* storage, uint32_t numSamplesInPopulation){
	genometools::THardyWeinbergGenotypeProbabilities pGenotype;
	for(int i=0; i<initialGridSize; ++i){
		pGenotype.set(_initialGrid[i].f);
		_initialGrid[i].density = _calcPosterior(storage, numSamplesInPopulation, pGenotype);
	}
};

void TAlleleFreqEstimatorBayes::_estimateMAP(const TSampleLikelihoods* storage, uint32_t numSamplesInPopulation){
	//use simple line-search to find MAP
	//initialize
	genometools::THardyWeinbergGenotypeProbabilities pGenotype;
	pGenotype.set(_guessInitialAlleleFrequency(storage, numSamplesInPopulation));
	logDensity_atMAP = _calcPosterior(storage, numSamplesInPopulation, pGenotype);
	double step = 0.05;
	if(pGenotype.f() > 0.5){
		step = -0.05;
	}

	//now do line search
	for(int i=0; i<numMAPSIterations; ++i){
		//propose new f, ensure [0,1] interval
		if(pGenotype.f() + step < 0.0){
			pGenotype.set(0.0);
		} else if(pGenotype.f() + step > 1.0){
			pGenotype.set(1.0);
		} else {
			pGenotype.set(pGenotype.f() + step);
		}

		//calc LL and switch if LL is lower
		double LL = _calcPosterior(storage, numSamplesInPopulation, pGenotype);
		if(LL < logDensity_atMAP || pGenotype.f() == 0.0 || pGenotype.f() == 1.0){
			step = - step / 2.718282;
		}
		logDensity_atMAP = LL;
	}
	f_MAP = pGenotype.f();

	//assumes that initial grid was pre-calculated! to check and boundaries
	_fillInitialGrid(storage, numSamplesInPopulation);

	//check if MAP is zero or one
	if(_initialGrid[0].density >= logDensity_atMAP){
		f_MAP = 0.0;
		logDensity_atMAP = _initialGrid[0].density;
	} else if(_initialGrid[initialGridLast].density >= logDensity_atMAP){
		f_MAP = 1.0;
		logDensity_atMAP = _initialGrid[initialGridLast].density;
	}
};

void TAlleleFreqEstimatorBayes::_estimateCredibleIntervals(const TSampleLikelihoods* storage, uint32_t numSamplesInPopulation){
	//use initial grid to define final grid
	//search first that is larger than LL_atMAP - logGridThreshold
	int first = 0;
	double relevantLL = logDensity_atMAP - logGridThreshold;
	while(_initialGrid[first].density < relevantLL && _initialGrid[first].f < f_MAP){
		++first;
	}
	if(first > 0) --first;

	//search last
	int last = initialGridLast;
	while(_initialGrid[last].density < relevantLL && _initialGrid[last].f > f_MAP){
		--last;
	}
	if(last < initialGridLast) ++last;

	//now prepare grid and calculate LL
	//make sure none is > LL_atMAP
	double step = (_initialGrid[last].f - _initialGrid[first].f) / (double) (gridLast);
	double integral = 0.0;
	genometools::THardyWeinbergGenotypeProbabilities pGenotype;
	for(int i=0; i<gridSize; ++i){
		_grid[i].f = _initialGrid[first].f + i * step;
		pGenotype.set(_grid[i].f);
		_grid[i].density = exp(_calcPosterior(storage, numSamplesInPopulation, pGenotype) - logDensity_atMAP);
		integral += _grid[i].density;
	}

	//adjust integral: remove half of first and last and multiply by step
	integral -= (_grid[0].density + _grid[gridLast].density) / 2.0; //first and last count half a step
	integral *= step;

	//find index left and right of MAP
	int left = gridLast;
	while(_grid[left].f >= f_MAP && left > 0){
		--left;
	}

	int right = 0;
	while(_grid[right].f <= f_MAP && right < gridLast){
		++right;
	}

	//add part from MAP to left and right grid points: average height as used when computing integral
	//then move left and right to next ones
	double halfStep = step / 2.0;
	double CI = 0;
	if(left > 0){
		CI += _grid[left].density * halfStep;
		-- left;
	}
	if(right < gridLast){
		CI += _grid[right].density * halfStep;
		++right;
	}

	//now find 90% CI by iteratively adding on left and right of MAP, depending on which has higher LL
	double relevantIntegral = credibleInterval * integral;
	if(CI >= relevantIntegral){
		f_CI_upper = _grid[right].f;
		f_CI_lower = _grid[left].f;
	} else {
		while(CI < relevantIntegral){
			if(right < gridSize && (left < 0 || _grid[left].density < _grid[right].density )){
				//add at right
				double add = (_grid[right].density + _grid[right - 1].density) * halfStep;
				if(CI + add >= relevantIntegral){
					f_CI_lower = _grid[left + 1].f;
					f_CI_upper = _grid[right - 1].f + step * (relevantIntegral - CI) / add;
					break;
				} else {
					CI += add;
					++right;
				}
			} else {
				//add at left
				double add = (_grid[left].density + _grid[left + 1].density) * halfStep;
				if(CI + add >= relevantIntegral){
					f_CI_upper = _grid[right - 1].f;
					f_CI_lower = _grid[left + 1].f + step * (relevantIntegral - CI) / add;
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

Probability TAlleleFreqEstimatorBayes::estimate(const TSampleLikelihoods* storage, uint32_t numSamplesInPop){
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

void TAlleleFreqEstimatorBayes::estimateAndWrite(const TSampleLikelihoods* storage, uint32_t numSamplesInPop, TOutputFile & out){
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
	genometools::THardyWeinbergGenotypeProbabilities pGenotype(f_MAP);
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
void TAlleleFreqMCMCOutput::initialize(std::string popString, genometools::TPopulationSamples & samples, std::string OutputName, TLog* logfile){
	//parse string to identify pops for which MCMC shoudl be written
	std::vector<std::string> tmp;
	coretools::str::fillContainerFromString(popString, tmp, ',');
	outputName = OutputName;

	//extract indexes
	if(tmp.size() > 0){
		//write all?
		if(tmp.size() == 1 && tmp[0] == "all"){
			for(size_t p=0; p<samples.numPopulations(); p++){
				popIndex.push_back(p);
				header.push_back(samples.getPopulationName(p));
			}
		} else {
			for(auto& name : tmp){
				if(samples.populationExists(name)){
					popIndex.push_back(samples.populationIndex(name));
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
TAlleleFreqEstimator::TAlleleFreqEstimator(TParameters &, TLog* Logfile){
	vcfRead = false;
	logfile = Logfile;
};

std::vector<std::string> TAlleleFreqEstimator::_composeHeaderAlleleFreq(bool writeGenoFreq, bool doBayesian, TAlleleFreqEstimatorBayes* BHWEstimator){
	std::vector<std::string> header = {"chr", "pos", "ref", "alt"};

	for(size_t p=0; p<samples.numPopulations(); p++){
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

void TAlleleFreqEstimator::_writeBayesianEstimatesOnePop(TOutputFile & out, TSampleLikelihoods* theseSamples, uint32_t numSamples, TAlleleFreqEstimatorBayes* BHWEstimator){
	out << BHWEstimator->estimate(theseSamples, numSamples); //MAP
	out << BHWEstimator->lowerCredibleInterval();
	out << BHWEstimator->upperCredibleInterval();
};

void TAlleleFreqEstimator::_writeEstimatesOnePop(TOutputFile & out, genometools::TGenotypeFrequencies & genoFrequencies, double alleleFrequency, TSampleLikelihoods* theseSamples, uint32_t numSamples, TAlleleFreqEstimatorHardyWeinberg & MLHWEstimator, TAlleleFreqEstimatorBayes* BHWEstimator, double epsF, bool writeGenoFreq, bool doBayesian){
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
	out << MLHWEstimator.estimate(theseSamples, numSamples, epsF);

	//Bayesian estimation
	if(doBayesian){
		_writeBayesianEstimatesOnePop(out, theseSamples, numSamples, BHWEstimator);
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
	if (samples.hasSamples())
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
	TAlleleFreqEstimatorBayes* BHWEstimator = nullptr;
	bool doBayesian = false;
	if(Parameters.parameterExists("doBayesian")){
		doBayesian = true;
		BHWEstimator = new TAlleleFreqEstimatorBayes(Parameters, logfile, randomGenerator);
	}

	bool writeGenoFreq = Parameters.parameterExists("writeGenoFreq");

	//output file
	std::string tmp = coretools::str::extractBeforeLast(vcfFilename, ".vcf");
	std::string outputName = Parameters.getParameterWithDefault<std::string>("out", tmp) + "_alleleFreq.txt.gz";
	logfile->list("Will write allele frequencies to file '" + outputName + "'.");
	TOutputFile out(outputName);

	//write header
	out.writeHeader(_composeHeaderAlleleFreq(writeGenoFreq, doBayesian, BHWEstimator));

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    while(reader.readDataFromVCF(storage, samples)){
    	//print SNP
 		reader.writePosition(out);

 		//write estimates based on genoFrequencies (if only 1 pop, use the one of reader)
 		if(samples.numPopulations() == 1){
 			_writeEstimatesOnePop(out, *(reader.genotypeFrequencies()), reader.alleleFrequency(), storage.samples(), samples.numSamples(), MLHWEstimator, BHWEstimator, epsF, writeGenoFreq, doBayesian);
 		} else {
 			genometools::TGenotypeFrequencies genoFrequencies;
 	 		for(size_t p=0; p<samples.numPopulations(); p++){
 	 			genoFrequencies.estimate(&storage[samples.startIndex(p)], samples.numSamplesInPop(p), epsF);
 	 			_writeEstimatesOnePop(out, genoFrequencies, genoFrequencies.alleleFrequency(), &storage[samples.startIndex(p)], samples.numSamplesInPop(p), MLHWEstimator, BHWEstimator, epsF, writeGenoFreq, doBayesian);
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

	for(size_t p=0; p<samples.numPopulations(); p++){
		std::string popName = samples.getPopulationName(p);
		header.push_back("numDiploid_" + popName);
		header.push_back("numHaploid_" + popName);

		BHWEstimator.composeHeader(header, popName);

		header.push_back("MCMCAcceptanceRate_" + popName);
	}

    //pairwise comparisons
	for(size_t p1=0; p1<(samples.numPopulations()-1); ++p1){
		for(size_t p2 = p1+1; p2 < samples.numPopulations(); ++p2){
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
	TAlleleFreqEstimatorBayes BHWEstimator(Parameters, logfile, randomGenerator);

	//genotype frequencies estimator
	genometools::TGenotypeFrequencies genoFrequencies;

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
	for(size_t p=0; p<samples.numPopulations(); p++){
		mcmcChains[p].resize(numIterations);
	}

	//output file
	std::string tmp = coretools::str::extractBeforeLast(vcfFilename, ".vcf");
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
    while(reader.readDataFromVCF(storage, samples)){
    	//print SNP
 		reader.writePosition(out);

 		//run MCMC
 		logfile->listFlush("Running estimates for " + reader.chr() + ":" + toString(reader.position()) + " ...");
 		for(size_t p=0; p<samples.numPopulations(); p++){
			//write num samples with data
			genoFrequencies.estimate(&storage[samples.startIndex(p)], samples.numSamplesInPop(p), Parameters.getParameterWithDefault("epsF", 0.0000001));
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
 		for(size_t p1=0; p1<(samples.numPopulations()-1); ++p1){
			for(size_t p2 = p1+1; p2 < samples.numPopulations(); ++p2){
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

void TAlleleFreqEstimator::writeAlleleFrequencyLikelihoods(TParameters & Parameters, TRandomGenerator*){
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
	std::string tmp = coretools::str::extractBeforeLast(vcfFilename, ".vcf");
	std::string outputName = Parameters.getParameterWithDefault<std::string>("out", tmp) + "_alleleFreqLikelihoods";
	logfile->list("Will write allele frequencies to files '" + outputName + "[POP].txt.gz'.");

	std::vector<TOutputFile> out(samples.numPopulations());
	if(samples.numPopulations() == 1){
		out[0].open(outputName + ".txt.gz");
		out[0].writeHeader(header);
		logfile->list("Will write allele frequency likelihoods to file '" + out[0].name() + "'.");
	} else {
		logfile->startIndent("Will write allele frequency likelihoods to files:");
		for(size_t p=0; p<samples.numPopulations(); p++){
			out[p].open(outputName + "_" + samples.getPopulationName(p) + ".txt.gz");
			out[p].writeHeader(header);
			logfile->list(out[p].name());
		}
		logfile->endIndent();
	}

	//prepare genotype probability object
	genometools::THardyWeinbergGenotypeProbabilities genoProb;

	//run through VCF file
	logfile->startIndent("Parsing VCF file:");
	while(reader.readDataFromVCF(storage, samples)){
		//calculate and write allele frequency likelihoods for every population
		for(size_t p=0; p<samples.numPopulations(); p++){
			out[p] << reader.chr() << reader.position();
			for(auto& f : freq){
				genoProb.set(f);
				coretools::LogProbability LL = 0.0;
				for(uint32_t i=0; i<samples.numSamplesInPop(p); i++){
					LL += storage[samples.startIndex(p) + i].HWESum<coretools::LogProbability>(genoProb);
				}
				out[p] << LL;
			}
			out[p] << std::endl;
		}
	}

	//clean up
	for(size_t p=0; p<samples.numPopulations(); p++){
		out[p].close();
	}
	_closeVCF();

};

}; //end namespace
