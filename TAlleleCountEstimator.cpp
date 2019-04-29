/*
 * TAlleleCountEstimator.cpp
 *
 *  Created on: Dec 9, 2018
 *      Author: phaentu
 */


#include "TAlleleCountEstimator.h"
//-------------------------------------------------
// TSAFChooseStorage
//-------------------------------------------------
TSAFChooseStorage::TSAFChooseStorage(int k){
	_k = k;
	log_choose_k_j = new double[_k+1]; //from zero to k

	//now fill
	for(int j=0; j<=_k; j++)
		log_choose_k_j[j] = chooseLog(_k, j);
};

TSAFChooseStorage::~TSAFChooseStorage(){
	delete[] log_choose_k_j;
};

double TSAFChooseStorage::logChoose(int j){
	if(j > _k) throw "j must be <= k in TSAFChooseStorage!";
	return log_choose_k_j[j];
};

double TSAFChooseStorage::operator[](int j){
	return log_choose_k_j[j];
};

//-------------------------------------------------
// TSiteAlleleFrequencyLikelihoods
//-------------------------------------------------
TSiteAlleleFrequencyLikelihoods::TSiteAlleleFrequencyLikelihoods(int numIndividuals){
	storageSize = 0;
	numAlleleCounts = 2*numIndividuals + 1;
	updateStorage(numAlleleCounts);
	logOf2 = log(2.0);
};

TSiteAlleleFrequencyLikelihoods::~TSiteAlleleFrequencyLikelihoods(){
	delete[] log_alleleFrequencyLikelihoods_h;
};

void TSiteAlleleFrequencyLikelihoods::updateStorage(int numIndividuals){
	int numRequiredAlleleCounts = 2*numIndividuals + 1;
	if(numRequiredAlleleCounts > storageSize){
		delete[] log_alleleFrequencyLikelihoods_h;
		log_alleleFrequencyLikelihoods_h = new double[numRequiredAlleleCounts];
		storageSize = numRequiredAlleleCounts;
	}
};

TSAFChooseStorage* TSiteAlleleFrequencyLikelihoods::getLogChoose(int counts){
	//did we already calculate it?
	std::map<int, TSAFChooseStorage*>::iterator it = log_choose.find(counts);
	if(it == log_choose.end()){
		log_choose.emplace(counts, new TSAFChooseStorage(counts));
		it = log_choose.find(counts);
	}
	return it->second;
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

void TSiteAlleleFrequencyLikelihoods::fillLog(const TPopulationLikehoodSample* data, int numSamples){
	//Calculating allele frequency likelihoods according to Nielsen et al. (2012) PLoS One, page 3
	//adapted to also work for haploid individuals (which only have likelihoods for genotypes 0 and 1)
	//all calculations done in log

	//set all h_j = 0
	for(int j=0; j<numAlleleCounts; j++)
		log_alleleFrequencyLikelihoods_h[j] = 0.0;

	//initialize with first individual
	if(data[0].isHaploid){
		log_alleleFrequencyLikelihoods_h[0] = qualMap.phredIntToLogErrorMap[data[0].phredLikelihoods[0]];
		log_alleleFrequencyLikelihoods_h[1] = qualMap.phredIntToLogErrorMap[data[0].phredLikelihoods[1]];
		numAlleleCounts = 1;
	} else {
		log_alleleFrequencyLikelihoods_h[0] =          qualMap.phredIntToLogErrorMap[data[0].phredLikelihoods[0]];
		log_alleleFrequencyLikelihoods_h[1] = logOf2 + qualMap.phredIntToLogErrorMap[data[0].phredLikelihoods[1]];
		log_alleleFrequencyLikelihoods_h[2] =          qualMap.phredIntToLogErrorMap[data[0].phredLikelihoods[2]];
		numAlleleCounts = 2;
	}

	//Recursion
	for(int s=1; s<numSamples; ++s){
		if(!data[0].isMissing){
			int j = numAlleleCounts;

			if(data[0].isHaploid){
				//first fill new ones to avoid multiplication with zero (relevant in log)
				log_alleleFrequencyLikelihoods_h[j+1] = qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[1]] + log_alleleFrequencyLikelihoods_h[j];

				//now fill those already used
				for(; j>1; j--){
					log_alleleFrequencyLikelihoods_h[j] = protectedSumInLog(
							qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[1]] + log_alleleFrequencyLikelihoods_h[j-1],
							qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[0]] + log_alleleFrequencyLikelihoods_h[j]    );
				}

				//special case for j=1,0
				log_alleleFrequencyLikelihoods_h[1] = protectedSumInLog(
							qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[1]] + log_alleleFrequencyLikelihoods_h[0],
							qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[0]] + log_alleleFrequencyLikelihoods_h[1]  );
				log_alleleFrequencyLikelihoods_h[0] = qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[0]] + log_alleleFrequencyLikelihoods_h[0];

				//increase total number of haplotypes by one
				numAlleleCounts += 1;
			} else {
				//first fill new ones to avoid multiplication with zero (relevant in log)
				log_alleleFrequencyLikelihoods_h[j+2] = qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[2]] + log_alleleFrequencyLikelihoods_h[j];
				log_alleleFrequencyLikelihoods_h[j+1] = protectedSumInLog(
															qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[2]] + log_alleleFrequencyLikelihoods_h[j-1],
												   logOf2 + qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[1]] + log_alleleFrequencyLikelihoods_h[j]    );

				//now fill those already used
				for(; j>1; j--){
					log_alleleFrequencyLikelihoods_h[j] = protectedSumInLog(
									 qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[2]] + log_alleleFrequencyLikelihoods_h[j-2],
							logOf2 + qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[1]] + log_alleleFrequencyLikelihoods_h[j-1],
									 qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[0]] + log_alleleFrequencyLikelihoods_h[j]    );
				}

				//special case for j=1,0
				log_alleleFrequencyLikelihoods_h[1] = protectedSumInLog(
							logOf2 + qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[1]] + log_alleleFrequencyLikelihoods_h[0],
									 qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[0]] + log_alleleFrequencyLikelihoods_h[1]  );

				log_alleleFrequencyLikelihoods_h[0] = qualMap.phredIntToLogErrorMap[data[s].phredLikelihoods[0]] + log_alleleFrequencyLikelihoods_h[0];

				//increase total number of haplotypes by two
				numAlleleCounts += 2;
			}
		}
	}

	//Termination: add binomial coefficient
	TSAFChooseStorage* logChoose = getLogChoose(numAlleleCounts);
	for(int j=0; j<numAlleleCounts; j++)
		log_alleleFrequencyLikelihoods_h[j] -= logChoose->logChoose(j);

	//Normalization
	normalize();
};

void TSiteAlleleFrequencyLikelihoods::fillNatural(const TPopulationLikehoodSample* data, int numSamples){
	//Calculating allele frequency likelihoods according to Nielsen et al. (2012) PLoS One, page 3
	//adapted to also work for haploid individuals (which only have likelihoods for genotypes 0 and 1)

	//set all h_j = 0
	for(int j=0; j<numAlleleCounts; j++)
		log_alleleFrequencyLikelihoods_h[j] = 0.0;

	//initialize with first individual
	if(data[0].isHaploid){
		log_alleleFrequencyLikelihoods_h[0] =     qualMap.phredIntToErrorMap[data[0].phredLikelihoods[0]];
		log_alleleFrequencyLikelihoods_h[1] =     qualMap.phredIntToErrorMap[data[0].phredLikelihoods[1]];
		numAlleleCounts = 1;
	} else {
		log_alleleFrequencyLikelihoods_h[0] =     qualMap.phredIntToErrorMap[data[0].phredLikelihoods[0]];
		log_alleleFrequencyLikelihoods_h[1] = 2 * qualMap.phredIntToErrorMap[data[0].phredLikelihoods[1]];
		log_alleleFrequencyLikelihoods_h[2] =     qualMap.phredIntToErrorMap[data[0].phredLikelihoods[2]];
		numAlleleCounts = 2;
	}

	//Recursion
	for(int s=1; s<numSamples; ++s){
		if(!data[0].isMissing){
			int j = numAlleleCounts;

			if(data[0].isHaploid){
				//first fill new ones to avoid multiplication with zero (relevant in log, kept here to code consistent)
				log_alleleFrequencyLikelihoods_h[j+1] = qualMap.phredIntToErrorMap[data[s].phredLikelihoods[2]] * log_alleleFrequencyLikelihoods_h[j-1];

				//now fill those already used
				for(; j>1; j--){
					log_alleleFrequencyLikelihoods_h[j] = qualMap.phredIntToErrorMap[data[s].phredLikelihoods[1]] * log_alleleFrequencyLikelihoods_h[j-1]
														+ qualMap.phredIntToErrorMap[data[s].phredLikelihoods[0]] * log_alleleFrequencyLikelihoods_h[j];
				}

				//special case for j=1,0
				log_alleleFrequencyLikelihoods_h[1] = qualMap.phredIntToErrorMap[data[s].phredLikelihoods[1]] * log_alleleFrequencyLikelihoods_h[0]
													+ qualMap.phredIntToErrorMap[data[s].phredLikelihoods[0]] * log_alleleFrequencyLikelihoods_h[1];
				log_alleleFrequencyLikelihoods_h[0] = qualMap.phredIntToErrorMap[data[s].phredLikelihoods[0]] * log_alleleFrequencyLikelihoods_h[0];

				//increase total number of haplotypes by one
				numAlleleCounts += 1;
			} else {
				//first fill new ones to avoid multiplication with zero (relevant in log, kept here to code consistent)
				log_alleleFrequencyLikelihoods_h[j+2] =     qualMap.phredIntToErrorMap[data[s].phredLikelihoods[2]] * log_alleleFrequencyLikelihoods_h[j];
				log_alleleFrequencyLikelihoods_h[j+1] =     qualMap.phredIntToErrorMap[data[s].phredLikelihoods[2]] * log_alleleFrequencyLikelihoods_h[j-1]
													  + 2 * qualMap.phredIntToErrorMap[data[s].phredLikelihoods[1]] * log_alleleFrequencyLikelihoods_h[j];

				//now fill those already used
				for(; j>1; j--){
					log_alleleFrequencyLikelihoods_h[j] =     qualMap.phredIntToErrorMap[data[s].phredLikelihoods[2]] * log_alleleFrequencyLikelihoods_h[j-2]
														+ 2 * qualMap.phredIntToErrorMap[data[s].phredLikelihoods[1]] * log_alleleFrequencyLikelihoods_h[j-1]
														+     qualMap.phredIntToErrorMap[data[s].phredLikelihoods[0]] * log_alleleFrequencyLikelihoods_h[j];
				}

				//special case for j=1,0
				log_alleleFrequencyLikelihoods_h[1] = 2 * qualMap.phredIntToErrorMap[data[s].phredLikelihoods[1]] * log_alleleFrequencyLikelihoods_h[0]
													+     qualMap.phredIntToErrorMap[data[s].phredLikelihoods[0]] * log_alleleFrequencyLikelihoods_h[1];
				log_alleleFrequencyLikelihoods_h[0] =     qualMap.phredIntToErrorMap[data[s].phredLikelihoods[0]] * log_alleleFrequencyLikelihoods_h[0];

				//increase total number of haplotypes by two
				numAlleleCounts += 2;
			}
		}
	}

	//Termination: put in log and add binomial coefficient
	TSAFChooseStorage* logChoose = getLogChoose(numAlleleCounts);
	for(int j=0; j<numAlleleCounts; j++)
		log_alleleFrequencyLikelihoods_h[j] = log(log_alleleFrequencyLikelihoods_h[j]) - logChoose->logChoose(j);

	//Normalization
	normalize();
};

void TSiteAlleleFrequencyLikelihoods::fill(const TPopulationLikehoodSample* data, int numSamples){
	//smallest likelihood is 10^-25.5 (phred 255).
	//A double can store up to 10^-308.
	//Hence we can store up to (10^25.5)^12 without underflow
	if(numSamples > 12)
		fillLog(data, numSamples);
	else
		fillNatural(data, numSamples);
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
	TPopulationLikelihoodReader reader(params, logfile, false);
	reader.openVCF(vcfFilename, logfile);
	logfile->endIndent();

	//Match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		 samples.readSamplesFromVCFNames(reader.getSampleVCFNames());

	//prepare site allele frequency likelihood calculators
	TSiteAlleleFrequencyLikelihoods** saf = new TSiteAlleleFrequencyLikelihoods*[samples.numPopulations()];
	for(int p=0; p<samples.numPopulations(); p++){
		saf[p] = new TSiteAlleleFrequencyLikelihoods(samples.numSamplesInPop(p));
	}

	//open output file
	std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
	std::string outname = params.getParameterStringWithDefault("out", tmp);
	std::string filename = outname + "_alleleCounts.txt.gz";
	logfile->list("Will write estimated allele counts to file '" + outname + "'.");
	gz::ogzstream aleleCountFile(filename.c_str());
	if(!aleleCountFile)
		throw "Failed to open file '" + filename + "' for writing!";

	//write header
	bool useLocusName = params.parameterExists("useLocusName");
	char sep = '\t';
	if(useLocusName){
		logfile->list("Will print locus names (rather than chromosome and position).");
		aleleCountFile << "Locus";
		sep = '_';
	} else
		aleleCountFile << "chr\tpos";
	for(int p=0; p<samples.numPopulations(); p++)
		aleleCountFile << "\t" << samples.getPopulationName(p);
	aleleCountFile << "\n";

	// initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);
	TPopulationLikehoodStorage data(samples.numSamples());

	//run through VCF file
	logfile->startIndent("Parsing VCF file and estimating allele counts:");
	while(reader.readDataFromVCF(data, samples, logfile)){
		//write chromosome and position
		aleleCountFile << reader.chr() << sep << reader.position();

		//print MLE count for each population
		for(int p=0; p<samples.numPopulations(); p++){
			//calculate allele frequency likelihoods
			saf[p]->fill(&data[samples.startIndex(p)], samples.numSamplesInPop(p));

			//and print MLE counts
			//TODO: find way to estimate counts among samples with data!
			//aleleCountFile << "\t" << saf[p]->getMLAlleleCount(*randomGenerator) << "/" << 2*samples.numSamplesWithDataInPop(sampleIsMissing, p);
			aleleCountFile << "\t" << saf[p]->getMLAlleleCount(*randomGenerator) << "/" << 2 * samples.numSamplesInPop(p);
		}
		aleleCountFile << std::endl;
	}

	//clean up
	for(int p=0; p<samples.numPopulations(); p++)
		delete saf[p];
	delete[] saf;

	//report final status
	logfile->endIndent();
	reader.concludeFilters(logfile);
	logfile->endIndent();
};
