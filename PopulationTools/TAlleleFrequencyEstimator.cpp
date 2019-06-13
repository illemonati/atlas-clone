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
	double mcmcLength = Parameters.getParameterIntWithDefault("mcmcLength", 100000);
	logfile->list("Will run MCMC chains of length " + toString(mcmcLength) + ".");
	numBurnins = Parameters.getParameterIntWithDefault("numBurnins", 3);
	burninLength = Parameters.getParameterIntWithDefault("burninLength", 1000);
	logfile->list("Will run " + toString(numBurnins) + " burnins of length " + toString(burninLength) + " each.");

	alpha = Parameters.getParameterDoubleWithDefault("alpha", 0.5);
	beta = Parameters.getParameterDoubleWithDefault("beta", 0.5);
	oneMinusAlpha = 1.0 - alpha;
	oneMinusBeta = 1.0 - beta;
	logfile->list("Will use a beta(" + toString(alpha) + "," + toString(beta) + ") prior.");
	logfile->endIndent();

	//initialize
	mcmc.resize(mcmcLength);
	lowerIndex = floor(mcmcLength * 0.05);
	upperIndex = ceil(mcmcLength * 0.95);
	f = 0.5;
	f_lower = 0.0;
	f_upper = 1.0;
};

double TAlleleFreqEstimatorBayes::guessInitialAlleleFrequency(TPopulationLikehoodLocus & storage, TGlfConverter & glfConverter){
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

int TAlleleFreqEstimatorBayes::makeMCMCUpdate(TPopulationLikehoodLocus & storage, double & oldLL, const double & prop, THardyWeinbergGenotypeProbabilities* pGenotype, int & old, TGlfConverter & glfConverter){
	int cur = 1-old;

	//propose new f
	double new_f = fabs(pGenotype[old].f + randomGenerator->getRand() * prop - prop / 2.0);
	if(new_f > 1.0) new_f = 2.0 - new_f;
	pGenotype[cur].set(new_f);

	//calc hastings
	double LL = calcLL(storage, pGenotype[cur], glfConverter);
	double priorRatio = oneMinusAlpha * (log(pGenotype[cur].f) - log(pGenotype[old].f)) + oneMinusBeta * (log(pGenotype[old].oneMinusf) - log(pGenotype[old].oneMinusf));
	double hastings = LL - oldLL + priorRatio;

	//accept or reject
	if(log(randomGenerator->getRand()) < hastings){
		old = cur;
		oldLL = LL;
		return true;
	}
	return false;
};

double TAlleleFreqEstimatorBayes::estimate(TPopulationLikehoodLocus & storage, TGlfConverter & glfConverter){
	//initialize
	THardyWeinbergGenotypeProbabilities pGenotype[2];
	int old = 0;
	pGenotype[old].set(guessInitialAlleleFrequency(storage, glfConverter));
	double prop = pGenotype[old][1] * 0.1;

	//calc initial LL
	double oldLL = calcLL(storage, pGenotype[old], glfConverter);

	//run burnins
	for(int b=0; b<numBurnins; b++){
		int acceptance = 0;
		for(int i=0; i<burninLength; i++){
			acceptance += makeMCMCUpdate(storage, oldLL, prop, pGenotype, old, glfConverter);
		}

		//adjust proposal range
		prop = prop * 3.0 * (double) acceptance / (double) burninLength;
	}

	//run MCMC
	for(size_t i=0; i<mcmc.size(); i++){
		makeMCMCUpdate(storage, oldLL, prop, pGenotype, old, glfConverter);

		//store current f
		mcmc[i] = pGenotype[1-old].f;
	}

	//estimate posterior mean
	double f = 0.0;
	for(const double& sample : mcmc){
		f += sample;
	}
	f /= (double) mcmc.size();

	//estimate posterior quantiles
	std::sort(mcmc.begin(), mcmc.end());
	f_lower = mcmc[lowerIndex];
	f_upper = mcmc[upperIndex];

	return f;
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
	std::vector<std::string> header = {"chr", "pos", "ref", "alt", "numDiploid", "numHaploid", "freqAltHW", "freqGenoRefRef", "freqGenoRefAlt", "freqGenoAltAlt", "freqGenoRef", "freqGenoAlt", "freqAltGeno"};
	if(doBayesian){
		header.push_back("freqAltHWBayes");
		header.push_back("freqAltHWBayes_CI0.05");
		header.push_back("freqAltHWBayes_CI0.95");
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

    //report final status
	logfile->endIndent();
	reader.concludeFilters(logfile);
	if(reader.numAcceptedLoci() < 1)
		throw "No usable loci in VCF file '" + vcfFilename + "'!";
	logfile->endIndent();
};

