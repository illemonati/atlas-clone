/*
 * TAlleleCountEstimator.cpp
 *
 *  Created on: Dec 9, 2018
 *      Author: phaentu
 */


#include "TAlleleCountEstimator.h"

//-------------------------------------------------
// TSiteAlleleFrequencyLikelihoods
//-------------------------------------------------
TSiteAlleleFrequencyLikelihoods::TSiteAlleleFrequencyLikelihoods(int numIndividuals){
	numInd_k = numIndividuals;
	numAlleleCounts = 2*numInd_k + 1;
	log_alleleFrequencyLikelihoods_h = new double[numAlleleCounts];
	logOf2 = log(2.0);

	//fill inverse of choose
	log_choose_2k_j = new double[numAlleleCounts];
	for(int j=0; j<numAlleleCounts; ++j)
		log_choose_2k_j[j] = chooseLog(2*numInd_k, j);
};

TSiteAlleleFrequencyLikelihoods::~TSiteAlleleFrequencyLikelihoods(){
	delete[] log_alleleFrequencyLikelihoods_h;
	delete[] log_choose_2k_j;
};

double TSiteAlleleFrequencyLikelihoods::protectedSumInLog(double a, double b){
  //returns log(exp(a)+exp(b)) while protecting for underflow, inspired by ANGSD
  double maxVal;
  if(a>b) maxVal = a;
  else maxVal = b;
  double sumVal = exp(a-maxVal)+exp(b-maxVal);
  return log(sumVal) + maxVal;
};

double TSiteAlleleFrequencyLikelihoods::protectedSumInLog(double a, double b, double c){
  //returns log(exp(a)+exp(b)+exp(c)) while protecting for underflow, inspired by ANGSD
  double maxVal;
  if(a > b && a > c) maxVal = a;
  else if(b > c) maxVal = b;
  else maxVal = c;
  double sumVal = exp(a-maxVal)+exp(b-maxVal)+exp(c-maxVal);
  return log(sumVal) + maxVal;
};

void TSiteAlleleFrequencyLikelihoods::normalize(){
	double max = log_alleleFrequencyLikelihoods_h[0];
	for(int j=1; j<numAlleleCounts; j++){
		if(log_alleleFrequencyLikelihoods_h[j] > max)
			max = log_alleleFrequencyLikelihoods_h[j];
	}
	for(int j=0; j<numAlleleCounts; j++)
		log_alleleFrequencyLikelihoods_h[j] -= max;
};

void TSiteAlleleFrequencyLikelihoods::fillLog(uint8_t* phred){
	//Calculating allele frequency likelihoods according to Nielsen et al. (2012) PLoS One, page 3
	//initialize
	log_alleleFrequencyLikelihoods_h[0] = qualMap.phredIntToLogErrorMap[phred[0]];
	log_alleleFrequencyLikelihoods_h[1] = logOf2 + qualMap.phredIntToLogErrorMap[phred[1]];
	log_alleleFrequencyLikelihoods_h[2] = qualMap.phredIntToLogErrorMap[phred[2]];

	for(int j=3; j<numAlleleCounts; j++)
		log_alleleFrequencyLikelihoods_h[j] = 0.0;

	//Recursion
	for(int d=1; d<numInd_k; ++d){
		int s = 3*d;
		int j=2*d;

		//first fill new ones to avoid multiplication with zero (relevant in log)
		log_alleleFrequencyLikelihoods_h[j+2] = qualMap.phredIntToLogErrorMap[phred[s + 2]] + log_alleleFrequencyLikelihoods_h[j];
		log_alleleFrequencyLikelihoods_h[j+1] = protectedSumInLog(
													qualMap.phredIntToLogErrorMap[phred[s + 2]] + log_alleleFrequencyLikelihoods_h[j-1],
										   logOf2 + qualMap.phredIntToLogErrorMap[phred[s + 1]] + log_alleleFrequencyLikelihoods_h[j]    );

		//now fill those already used
		for(; j>1; j--){
			log_alleleFrequencyLikelihoods_h[j] = protectedSumInLog(
							 qualMap.phredIntToLogErrorMap[phred[s + 2]] + log_alleleFrequencyLikelihoods_h[j-2],
					logOf2 + qualMap.phredIntToLogErrorMap[phred[s + 1]] + log_alleleFrequencyLikelihoods_h[j-1],
							 qualMap.phredIntToLogErrorMap[phred[s]]     + log_alleleFrequencyLikelihoods_h[j]    );
		}

		//special case for j=1,0
		log_alleleFrequencyLikelihoods_h[1] = protectedSumInLog(
					logOf2 + qualMap.phredIntToLogErrorMap[phred[s + 1]] + log_alleleFrequencyLikelihoods_h[0],
							 qualMap.phredIntToLogErrorMap[phred[s]]     + log_alleleFrequencyLikelihoods_h[1]  );

		log_alleleFrequencyLikelihoods_h[0] = qualMap.phredIntToLogErrorMap[phred[s]] + log_alleleFrequencyLikelihoods_h[0];
	}

	//Termination
	for(int j=0; j<numAlleleCounts; j++)
		log_alleleFrequencyLikelihoods_h[j] -= log_choose_2k_j[j];

	//Normalization
	normalize();
};

void TSiteAlleleFrequencyLikelihoods::fillNatural(uint8_t* phred){
	//Calculating allele frequency likelihoods according to Nielsen et al. (2012) PLoS One, page 3
	//initialize
	log_alleleFrequencyLikelihoods_h[0] = qualMap.phredIntToErrorMap[phred[0]];
	log_alleleFrequencyLikelihoods_h[1] = 2 * qualMap.phredIntToErrorMap[phred[1]];
	log_alleleFrequencyLikelihoods_h[2] = qualMap.phredIntToErrorMap[phred[2]];

	for(int j=3; j<numAlleleCounts; j++)
		log_alleleFrequencyLikelihoods_h[j] = 0.0;

	//Recursion
	for(int d=1; d<numInd_k; ++d){
		int s = 3*d;
		int j=2*d;

		//first fill new ones to avoid multiplication with zero (relevant in log)
		log_alleleFrequencyLikelihoods_h[j+2] = qualMap.phredIntToErrorMap[phred[s + 2]] * log_alleleFrequencyLikelihoods_h[j];
		log_alleleFrequencyLikelihoods_h[j+1] = qualMap.phredIntToErrorMap[phred[s + 2]] * log_alleleFrequencyLikelihoods_h[j-1]
											  + 2 * qualMap.phredIntToErrorMap[phred[s + 1]] * log_alleleFrequencyLikelihoods_h[j];

		//now fill those already used
		for(; j>1; j--){
			log_alleleFrequencyLikelihoods_h[j] = qualMap.phredIntToErrorMap[phred[s + 2]] * log_alleleFrequencyLikelihoods_h[j-2]
												+ 2 * qualMap.phredIntToErrorMap[phred[s + 1]] * log_alleleFrequencyLikelihoods_h[j-1]
												+ qualMap.phredIntToErrorMap[phred[s]]     * log_alleleFrequencyLikelihoods_h[j];
		}

		//special case for j=1,0
		log_alleleFrequencyLikelihoods_h[1] = 2 * qualMap.phredIntToErrorMap[phred[s + 1]] * log_alleleFrequencyLikelihoods_h[0]
											+ qualMap.phredIntToErrorMap[phred[s]]     * log_alleleFrequencyLikelihoods_h[1];

		log_alleleFrequencyLikelihoods_h[0] = qualMap.phredIntToErrorMap[phred[s]] * log_alleleFrequencyLikelihoods_h[0];
	}

	//Termination
	for(int j=0; j<numAlleleCounts; j++)
		log_alleleFrequencyLikelihoods_h[j] = log(log_alleleFrequencyLikelihoods_h[j]) - log_choose_2k_j[j];

	//Normalization
	normalize();
};

void TSiteAlleleFrequencyLikelihoods::fill(uint8_t* phred){
	//smallest likelihood is 10^-25.5 (phred 255).
	//A double can store up to 10^-308.
	//Hence we can store up to (10^25.5)^12 without underflow
	if(numInd_k > 12)
		fillLog(phred);
	else
		fillNatural(phred);
};

int TSiteAlleleFrequencyLikelihoods::getMLAlleleCount(TRandomGenerator & randomGenerator){
	//first find ML and store all indexes that are at ML
	double ML = log_alleleFrequencyLikelihoods_h[0];
	for(int j=1; j<numAlleleCounts; j++){
		if(log_alleleFrequencyLikelihoods_h[j] > ML)
			ML = log_alleleFrequencyLikelihoods_h[j];
	}

	//now store all index at ML
	std::vector<int> MLEs;
	for(int j=0; j<numAlleleCounts; j++){
		if(log_alleleFrequencyLikelihoods_h[j] == ML){
			MLEs.emplace_back(j);
		}
	}

	//now choose randomly among those ate MLE
	return MLEs[randomGenerator.pickOne(MLEs.size())];
};

void TSiteAlleleFrequencyLikelihoods::print(){
	for(int j=0; j<numAlleleCounts; j++){
		std::cout << "\t" << log_alleleFrequencyLikelihoods_h[j];
	}
};


//-------------------------------------------------
// TAlleleCountEstimator
//-------------------------------------------------
TAlleleCountEstimator::TAlleleCountEstimator(TParameters & params, TLog* Logfile){
	logfile = Logfile;

	//initialize random generator
	//TODO: do the random generator initialization in the task switcher?
	logfile->listFlush("Initializing random generator ...");
	if(params.parameterExists("fixedSeed")){
		randomGenerator = new TRandomGenerator(params.getParameterLong("fixedSeed"), true);
	} else if(params.parameterExists("addToSeed")){
		randomGenerator = new TRandomGenerator(params.getParameterLong("addToSeed"), false);
	} else randomGenerator = new TRandomGenerator();
	logfile->write(" done with seed " + toString(randomGenerator->usedSeed) + "!");
};

TAlleleCountEstimator::~TAlleleCountEstimator(){
	delete randomGenerator;
};

void TAlleleCountEstimator::estimateAlleleCounts(TParameters & params){
	//read samples
	TPopulationSamples samples;
	if(params.parameterExists("samples"))
		samples.readSamples(params.getParameterString("samples"), logfile);

	//open VCF reader
	std::string vcfFilename = params.getParameterString("vcf");
	logfile->startIndent("Reading genotype likelihoods from VCF file '" + vcfFilename + "':");
	TPopulationLikelihoodReader reader(params, logfile);
	reader.openVCF(vcfFilename, logfile);
	logfile->endIndent();

	//Match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		 samples.readSamplesFromVCFNames(reader.getSampleVCFNames());

	// initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);
	uint8_t* curLocus = new uint8_t[samples.numSamples() * 3];

	//prepare site allele frequency likelihood calculators
	TSiteAlleleFrequencyLikelihoods** saf = new TSiteAlleleFrequencyLikelihoods*[samples.numPopulations()];
	for(int p=0; p<samples.numPopulations(); p++){
		saf[p] = new TSiteAlleleFrequencyLikelihoods(samples.numSamplesInPop(p));
	}

	//open output file
	std::string outname = params.getParameterStringWithDefault("out", "ATLAS_");
	std::string filename = outname + "alleleCounts.txt.gz";
	logfile->list("Will write estimated allele counts to file '" + outname + "'.");
	gz::ogzstream aleleCountFile(filename.c_str());
	if(!aleleCountFile)
		throw "Failed to open file '" + filename + "' for writing!";

	//write header
	aleleCountFile << "chr\tpos";
	for(int p=0; p<samples.numPopulations(); p++)
		aleleCountFile << "\t" << samples.getPopulationName(p);
	aleleCountFile << "\n";

	//run through VCF file
	logfile->startIndent("Parsing VCF file and estimating allele counts:");
	while(reader.readDataFromVCF(curLocus, samples, logfile)){
		//write chromosome and position
		aleleCountFile << reader.chr() << "\t" << reader.position();

		//print MLE count for each population
		for(int p=0; p<samples.numPopulations(); p++){
			//calculate allele frequency likelihoods
			saf[p]->fill(&curLocus[3*samples.startIndex(p)]);

			//and print MLE counts
			aleleCountFile << "\t" << saf[p]->getMLAlleleCount(*randomGenerator);
		}
		aleleCountFile << std::endl;
	}

	//clean up
	delete[] curLocus;
	for(int p=0; p<samples.numPopulations(); p++)
		delete saf[p];
	delete[] saf;

	//report final status
	logfile->endIndent();
	reader.concludeFilters(logfile);
	logfile->endIndent();
};
