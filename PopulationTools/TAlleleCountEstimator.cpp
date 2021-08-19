/*
 * TAlleleCountEstimator.cpp
 *
 *  Created on: Dec 9, 2018
 *      Author: phaentu
 */


#include "TAlleleCountEstimator.h"

namespace PopulationTools{

using genometools::homoFirst;
using genometools::het;
using genometools::homoSecond;
using genometools::haploidFirst;
using genometools::haploidSecond;

//-------------------------------------------------
// TSAFChooseStorage
//-------------------------------------------------
TSAFChooseStorage::TSAFChooseStorage(int k){
	_k = k;
	log_choose_k_j = new double[_k+1]; //from zero to k

	//now fill
	for(int j=0; j<=_k; j++){
		log_choose_k_j[j] = coretools::chooseLog(_k, j);
	}
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
	numAlleleCounts = 2*numIndividuals + 1;
	log_alleleFrequencyLikelihoods_h.resize(numAlleleCounts);
	logOf2 = log(2.0);
};

TSiteAlleleFrequencyLikelihoods::~TSiteAlleleFrequencyLikelihoods(){
	for(std::pair<int, TSAFChooseStorage*> it : log_choose){
		delete it.second;
	}
};

TSAFChooseStorage* TSiteAlleleFrequencyLikelihoods::_getLogChoose(int counts){
	//did we already calculate it?
	std::map<int, TSAFChooseStorage*>::iterator it = log_choose.find(counts);
	if(it == log_choose.end()){
		log_choose.emplace(counts, new TSAFChooseStorage(counts));
		it = log_choose.find(counts);
	}
	return it->second;
};

LogProbability TSiteAlleleFrequencyLikelihoods::_protectedSumInLog(const LogProbability a, const LogProbability b){
  //returns log(exp(a)+exp(b)) while protecting for underflow, inspired by ANGSD
  LogProbability maxVal;
  if(a>b) maxVal = a;
  else maxVal = b;
  Probability sumVal = (Probability) (a - maxVal) + (Probability) (b - maxVal);
  return LogProbability(sumVal) + maxVal;
};

LogProbability TSiteAlleleFrequencyLikelihoods::_protectedSumInLog(const LogProbability a, const LogProbability b, const LogProbability c){
  //returns log(exp(a)+exp(b)+exp(c)) while protecting for underflow, inspired by ANGSD
  LogProbability maxVal;
  if(a > b && a > c) maxVal = a;
  else if(b > c) maxVal = b;
  else maxVal = c;
  Probability sumVal = (Probability) (a - maxVal) + (Probability) (b - maxVal) + (Probability) (c - maxVal);
  return LogProbability(sumVal) + maxVal;
};

void TSiteAlleleFrequencyLikelihoods::normalize(){
	LogProbability max = log_alleleFrequencyLikelihoods_h[0];
	for(int j=1; j<numAlleleCounts+1; j++){
		if(log_alleleFrequencyLikelihoods_h[j] > max)
			max = log_alleleFrequencyLikelihoods_h[j];
	}
	for(int j=0; j<numAlleleCounts+1; j++)
		log_alleleFrequencyLikelihoods_h[j] -= max;
};

void TSiteAlleleFrequencyLikelihoods::_fillLog(TSampleLikelihoods* data, const uint32_t & numSamples){
	//Calculating allele frequency likelihoods according to Nielsen et al. (2012) PLoS One, page 3
	//adapted to also work for haploid individuals (which only have likelihoods for genotypes 0 and 1)
	//all calculations done in log

	//set all h_j = 0
	std::vector<double> alleleFrequencyLikelihoods_h(log_alleleFrequencyLikelihoods_h.size(), 0.0);
	numAlleleCounts = 0;

	//find first individual  with data
	int s = 0;
	while(s < numSamples && data[s].isMissing()){
		++s;
	}

	if(s < numSamples){
		//initialize with first individual
		if(data[s].isHaploid()){
			alleleFrequencyLikelihoods_h[0] = (LogProbability) data[s][haploidFirst];
			alleleFrequencyLikelihoods_h[1] = (LogProbability) data[s][haploidSecond];
			numAlleleCounts = 1;
		} else {
			alleleFrequencyLikelihoods_h[0] = (LogProbability) data[s][homoFirst];
			alleleFrequencyLikelihoods_h[1] = (LogProbability) data[s][het] + logOf2;
			alleleFrequencyLikelihoods_h[2] = (LogProbability) data[s][homoSecond];
			numAlleleCounts = 2;
		}

		//Recursion
		for(++s; s<numSamples; ++s){
			if(!data[s].isMissing()){
				int j = numAlleleCounts;

				if(data[s].isHaploid()){
					//first fill new ones to avoid multiplication with zero (relevant in log)
					alleleFrequencyLikelihoods_h[j+1] = (LogProbability) data[s][haploidSecond] + alleleFrequencyLikelihoods_h[j];

					//now fill those already used
					for(; j>0; j--){
						alleleFrequencyLikelihoods_h[j] = _protectedSumInLog(
								(LogProbability) data[s][haploidSecond] + alleleFrequencyLikelihoods_h[j-1],
								(LogProbability) data[s][haploidFirst] + alleleFrequencyLikelihoods_h[j]    );
					}

					//special case for j=0
					alleleFrequencyLikelihoods_h[0] = (LogProbability) data[s][haploidFirst] + alleleFrequencyLikelihoods_h[0];

					//increase total number of haplotypes by one
					numAlleleCounts += 1;
				} else {
					//first fill new ones to avoid multiplication with zero (relevant in log)
					alleleFrequencyLikelihoods_h[j+2] = (LogProbability) data[s][homoSecond] + alleleFrequencyLikelihoods_h[j];
					alleleFrequencyLikelihoods_h[j+1] = _protectedSumInLog(
																(LogProbability) data[s][homoSecond] + alleleFrequencyLikelihoods_h[j-1],
																(LogProbability) data[s][het] + alleleFrequencyLikelihoods_h[j] + logOf2 );

					//now fill those already used
					for(; j>1; j--){
						alleleFrequencyLikelihoods_h[j] = _protectedSumInLog(
										 (LogProbability) data[s][homoSecond] + alleleFrequencyLikelihoods_h[j-2],
										 (LogProbability) data[s][het] + alleleFrequencyLikelihoods_h[j-1] + logOf2,
										 (LogProbability) data[s][homoFirst] + alleleFrequencyLikelihoods_h[j] );
					}

					//special case for j=1,0
					alleleFrequencyLikelihoods_h[1] = _protectedSumInLog(
					        (LogProbability) data[s][het] + alleleFrequencyLikelihoods_h[0] + logOf2, // 2*het does not appear in equation in paper, but must be there
										 (LogProbability) data[s][homoFirst] + alleleFrequencyLikelihoods_h[1] );

					alleleFrequencyLikelihoods_h[0] = (LogProbability) data[s][homoFirst] + alleleFrequencyLikelihoods_h[0];

					//increase total number of haplotypes by two
					numAlleleCounts += 2;
				}
			}
		}

		//Termination: add binomial coefficient
		TSAFChooseStorage* logChoose = _getLogChoose(numAlleleCounts);
		for(int j=0; j<numAlleleCounts+1; j++){
			//numerical accuracy may raely lead to a value very slightly above 0.0. std::min is used to avoid an error when storing as LogProbability.
			log_alleleFrequencyLikelihoods_h[j] = std::min(alleleFrequencyLikelihoods_h[j] - logChoose->logChoose(j), 0.0);
		}

		//Normalization
		normalize();
	}
};

void TSiteAlleleFrequencyLikelihoods::_fillNatural(TSampleLikelihoods* data, const uint32_t & numSamples){
	//Calculating allele frequency likelihoods according to Nielsen et al. (2012) PLoS One, page 3
	//adapted to also work for haploid individuals (which only have likelihoods for genotypes 0 and 1)
	std::vector<double> alleleFrequencyLikelihoods_h(log_alleleFrequencyLikelihoods_h.size(), 0.0);

	numAlleleCounts = 0;

	//find first individual  with data
	int s = 0;
	while(s < numSamples && data[s].isMissing()){
		++s;
	}

	if(s < numSamples){
		//initialize with first individual
		if(data[s].isHaploid()){
			alleleFrequencyLikelihoods_h[0] = (Probability) data[s][haploidFirst];
			alleleFrequencyLikelihoods_h[1] = (Probability) data[s][haploidSecond];
			numAlleleCounts = 1;
		} else {
			alleleFrequencyLikelihoods_h[0] = (Probability) data[s][homoFirst];
			alleleFrequencyLikelihoods_h[1] = (Probability) data[s][het] * 2.0;
			alleleFrequencyLikelihoods_h[2] = (Probability) data[s][homoSecond];
			numAlleleCounts = 2;
		}

		//Recursion
		for(++s; s<numSamples; ++s){
			if(!data[s].isMissing()){
				int j = numAlleleCounts;

				if(data[s].isHaploid()){
					//first fill new ones to avoid multiplication with zero (relevant in log, kept here for code consistency)
					alleleFrequencyLikelihoods_h[j+1] = (Probability) data[s][haploidSecond] * alleleFrequencyLikelihoods_h[j];

					//now fill those already used
					for(; j>0; j--){
						alleleFrequencyLikelihoods_h[j] = (Probability) data[s][haploidSecond] * alleleFrequencyLikelihoods_h[j-1]
													    + (Probability) data[s][haploidFirst] * alleleFrequencyLikelihoods_h[j];
					}

					//special case for j=0
					alleleFrequencyLikelihoods_h[0] = (Probability) data[s][haploidFirst] * alleleFrequencyLikelihoods_h[0];

					//increase total number of haplotypes by one
					numAlleleCounts += 1;
				} else {
					//first fill new ones to avoid multiplication with zero (relevant in log, kept here to code consistent)
					alleleFrequencyLikelihoods_h[j+2] = (Probability) data[s][homoSecond] * alleleFrequencyLikelihoods_h[j];
					alleleFrequencyLikelihoods_h[j+1] = (Probability) data[s][homoSecond] * alleleFrequencyLikelihoods_h[j-1]
													  + 2.0 * (Probability) data[s][het] * alleleFrequencyLikelihoods_h[j];

					//now fill those already used
					for(; j>1; j--){
						alleleFrequencyLikelihoods_h[j] = (Probability) data[s][homoSecond] * alleleFrequencyLikelihoods_h[j-2]
														+ 2.0 * (Probability) data[s][het] * alleleFrequencyLikelihoods_h[j-1]
														+ (Probability) data[s][homoFirst] * alleleFrequencyLikelihoods_h[j];
					}

					//special case for j=1,0
					alleleFrequencyLikelihoods_h[1] = 2.0 * (Probability) data[s][het] * alleleFrequencyLikelihoods_h[0] // 2*het does not appear in equation in paper, but must be there
													+ (Probability) data[s][homoFirst] * alleleFrequencyLikelihoods_h[1];
					alleleFrequencyLikelihoods_h[0] = (Probability) data[s][homoFirst] * alleleFrequencyLikelihoods_h[0];

					//increase total number of haplotypes by two
					numAlleleCounts += 2;
				}
			}
		}

		//Termination: put in log and add binomial coefficient
		TSAFChooseStorage* logChoose = _getLogChoose(numAlleleCounts);
		for(int j=0; j<numAlleleCounts+1; j++){
			//numerical accuracy may raely lead to a value very slightly above 0.0. std::min is used to avoid an error when storing as LogProbability.
			log_alleleFrequencyLikelihoods_h[j] = std::min(log(alleleFrequencyLikelihoods_h[j]) - logChoose->logChoose(j), 0.0);
		}

		//Normalization
		normalize();
	}
};

void TSiteAlleleFrequencyLikelihoods::fill(TSampleLikelihoods* data, const uint32_t & numSamples){
	//smallest likelihood is 10^-25.5 (phred 255).
	//A double can store up to 10^-308.
	//Hence we can store up to (10^25.5)^12 without underflow
	if(numSamples > 12){
		_fillLog(data, numSamples);
	} else {
		_fillNatural(data, numSamples);
	}
};

void TSiteAlleleFrequencyLikelihoods::print(){
	for(int j=0; j<numAlleleCounts+1; j++){
		std::cout << "\t" << log_alleleFrequencyLikelihoods_h[j];
	}
};

void TSiteAlleleFrequencyLikelihoods::write(gz::ogzstream & file){
	for(int j=0; j<numAlleleCounts+1; j++){
		file << "\t" << log_alleleFrequencyLikelihoods_h[j];
	}
};

int TSiteAlleleFrequencyLikelihoods::getMLAlleleCount(coretools::TRandomGenerator & randomGenerator){
	//return 0 in case of no data
	if(numAlleleCounts == 0) return 0;

	//first find ML and store all indexes that are at ML
	double ML = log_alleleFrequencyLikelihoods_h[0];
	for(int j=1; j<numAlleleCounts+1; j++){
		if(log_alleleFrequencyLikelihoods_h[j] > ML)
			ML = log_alleleFrequencyLikelihoods_h[j];
	}

	//now store all index at ML
	std::vector<int> MLEs;
	for(int j=0; j<numAlleleCounts+1; j++){
		if(log_alleleFrequencyLikelihoods_h[j] == ML){
			MLEs.emplace_back(j);
		}
	}

	//now choose randomly among those ate MLE
	return MLEs[randomGenerator.sample(MLEs.size())];
};

const std::vector<coretools::LogProbability> & TSiteAlleleFrequencyLikelihoods::getLogAlleleFrequencyLikelihoods() const{
    return log_alleleFrequencyLikelihoods_h;
};

//-------------------------------------------------
// TAlleleCountEstimator
//-------------------------------------------------
TAlleleCountEstimator::TAlleleCountEstimator(coretools::TParameters & params, coretools::TLog* Logfile){
	logfile = Logfile;
};

TAlleleCountEstimator::~TAlleleCountEstimator(){};

void TAlleleCountEstimator::estimateAlleleCounts(coretools::TParameters & params, coretools::TRandomGenerator* randomGenerator){
	//read samples
	TPopulationSamples samples;
	if(params.parameterExists("samples"))
		samples.readSamples(params.getParameter<std::string>("samples"), logfile);

	//open VCF reader
	std::string vcfFilename = params.getParameter<std::string>("vcf");
	logfile->startIndent("Reading genotype likelihoods from VCF file '" + vcfFilename + "':");
	TPopulationLikelihoodReaderLocus reader(params, logfile, false);
	reader.openVCF(vcfFilename);
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

	//create out file
	std::string tmp = coretools::str::extractBeforeLast(vcfFilename, ".vcf");
	std::string outname = params.getParameterWithDefault<std::string>("out", tmp);
	std::string type = params.getParameterWithDefault<std::string>("outFormat", "default");
	TAlleleCountFile* alleleCountFile = prepareOutputFile(type, outname, params);

	//write header
	alleleCountFile->writeHeader(samples, params, logfile);


	// initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);
	TPopulationLikehoodLocus data(samples.numSamples());

	//run through VCF file
	logfile->startIndent("Parsing VCF file and estimating allele counts:");
	while(reader.readDataFromVCF(data, samples)){
		//write chromosome and position
		alleleCountFile->writePosition(reader.chr(), reader.position());

		//print MLE count for each population
		for(int p=0; p<samples.numPopulations(); p++){
			//calculate allele frequency likelihoods
			saf[p]->fill(&data[samples.startIndex(p)], samples.numSamplesInPop(p));

			//and print MLE counts
			alleleCountFile->writeCounts(saf[p]->getMLAlleleCount(*randomGenerator), saf[p]->getNumAlleles(), p);
		}
		alleleCountFile->endl();
	}

	//clean up
	for(int p=0; p<samples.numPopulations(); p++)
		delete saf[p];
	delete[] saf;
	delete alleleCountFile;

	//report final status
	logfile->endIndent();
	reader.concludeFilters();
	logfile->endIndent();
};

void TAlleleCountEstimator::writeAlleleFrequencyLikelihoods(coretools::TParameters & params){
	//TODO: write proper saf
	//read samples
	TPopulationSamples samples;
	if(params.parameterExists("samples"))
		samples.readSamples(params.getParameter<std::string>("samples"), logfile);

	//open VCF reader
	std::string vcfFilename = params.getParameter<std::string>("vcf");
	logfile->startIndent("Reading genotype likelihoods from VCF file '" + vcfFilename + "':");
	TPopulationLikelihoodReaderLocus reader(params, logfile, false);
	reader.openVCF(vcfFilename);
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
	std::string tmp = coretools::str::extractBeforeLast(vcfFilename, ".vcf");
	std::string outname = params.getParameterWithDefault<std::string>("out", tmp);
	std::string filename = outname + "_alleleFrequencyLikelihoods.txt.gz";
	logfile->list("Will write estimated allele counts to file '" + outname + "'.");
	gz::ogzstream alleleFrequencyLikelihoodFile(filename.c_str());
	if(!alleleFrequencyLikelihoodFile)
		throw "Failed to open file '" + filename + "' for writing!";

	//write header
	bool useLocusName = params.parameterExists("useLocusName");
	char sep = '\t';
	if(useLocusName){
		logfile->list("Will print locus names (rather than chromosome and position).");
		alleleFrequencyLikelihoodFile << "Locus";
		sep = '_';
	}
	alleleFrequencyLikelihoodFile << "chr" << sep << "pos\tnumAlleles";
	for(int p=0; p<samples.numPopulations(); p++)
		alleleFrequencyLikelihoodFile << "\t" << samples.getPopulationName(p);
	alleleFrequencyLikelihoodFile << "\n";

	// initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);
	TPopulationLikehoodLocus data(samples.numSamples());

	//run through VCF file
	logfile->startIndent("Parsing VCF file and estimating allele counts:");
	while(reader.readDataFromVCF(data, samples)){
		//write chromosome and position
		alleleFrequencyLikelihoodFile << reader.chr() << sep << reader.position();

		//print MLE count for each population
		for(int p=0; p<samples.numPopulations(); p++){
			//calculate allele frequency likelihoods
			saf[p]->fill(&data[samples.startIndex(p)], samples.numSamplesInPop(p));

			//print num alleles
			alleleFrequencyLikelihoodFile << "\t" << saf[p]->getNumAlleles();

			//print AFL
			saf[p]->write(alleleFrequencyLikelihoodFile);
		}
		alleleFrequencyLikelihoodFile << std::endl;
	}

	//clean up
	for(int p=0; p<samples.numPopulations(); p++)
		delete saf[p];
	delete[] saf;

	//report final status
	logfile->endIndent();
	reader.concludeFilters();
	logfile->endIndent();
};

TAlleleCountFile* TAlleleCountEstimator::prepareOutputFile(std::string type, std::string filePrefix, coretools::TParameters& params){
	//create output file
	TAlleleCountFile* alleleCountFile;
	if(type == "default"){
		filePrefix = filePrefix + "_alleleCounts.txt.gz";
		alleleCountFile = new TAlleleCountFile(filePrefix);
	} else if(type == "treemix"){
		filePrefix = filePrefix + "_treemix_alleleCounts.txt.gz";
		alleleCountFile = new TTreeMixFile(filePrefix);
	} else if(type == "flink"){
		filePrefix = filePrefix + "_flink_alleleCounts.txt.gz";
		alleleCountFile = new TFlinkFile(filePrefix);
	} else
		throw "Unknown output file type '" + type + "'!";

	logfile->list("Will write estimated allele counts to file '" + filePrefix + "'.");
	alleleCountFile->openFileToWrite(filePrefix);

	return(alleleCountFile);
}

void TAlleleCountEstimator::transformFormat(coretools::TParameters & params){
	// initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);

	//get parameters for in and output
	std::string countsFileName = params.getParameter<std::string>("alleleCounts");
	std::string tmp = coretools::str::extractBeforeLast(countsFileName, "_alleleCounts.txt.gz");
	std::string outname = params.getParameterWithDefault<std::string>("out", tmp);
	std::string type = params.getParameterWithDefault<std::string>("outFormat", "default");
	if(type == "default"){
		throw "Cannot transform alleleCounts file to original format!";
	}

	//open input file
	std::string origFileName = tmp + "_alleleCounts.txt.gz";
	gz::igzstream file(origFileName.c_str());
	if(!file) throw "Failed to open file '" + origFileName + " for reading!";

	//read header
	std::string line;
	std::getline(file, line);
	std::vector<std::string> tmp_vec;
	std::vector<std::string> populationNames;
	coretools::str::fillContainerFromStringWhiteSpace(line, tmp_vec, true);
	for(unsigned int i=2; i<tmp_vec.size(); ++i){
		populationNames.push_back(tmp_vec[i]);
	}

	//create output file
	TAlleleCountFile* alleleCountFile = prepareOutputFile(type, outname, params);

	//write header
	alleleCountFile->writeHeader(populationNames, params, logfile);

	//run through VCF file
	logfile->startIndent("Converting allele counts file:");
	while(file.good() && !file.eof()){
		std::vector<std::string> vec;
		std::getline(file, line);
		coretools::str::fillContainerFromStringWhiteSpace(line, vec, true);
		//write chromosome and position
		if(vec.size() > 0){
			alleleCountFile->writePosition(vec[0], vec[1]);

			//print MLE count for each population
			for(unsigned int p=2; p<vec.size(); p++){
				//print MLE counts
				std::vector<std::string> counts;
				coretools::str::fillContainerFromStringAny(vec[p], counts, "/");
				alleleCountFile->writeCounts(counts[0], counts[1], p-2);
			}
			alleleCountFile->endl();
		}
	}

	//close file
	file.close();

	//clean up
	delete alleleCountFile;

	//report final status
	logfile->endIndent();
};

}; //end namespace
