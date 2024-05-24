/*
 * TAlleleCountEstimator.cpp
 *
 *  Created on: Dec 9, 2018
 *      Author: phaentu
 */

#include "TAlleleCountEstimator.h"

#include "TAlleleCountFileFormat.h"
#include "TAlleleCountReader.h"
#include "TSafFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Math/mathFunctions.h"

namespace PopulationTools {

using coretools::LogProbability;
using coretools::Probability;
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::instances::randomGenerator;

namespace impl {

// Can't use Probability type as the sum inside the log will be > 1!
// Also, this function can return a value > 0 which will be normalized at the end
double _protectedSumInLog(double a, double b) {
	// returns log(exp(a)+exp(b)) while protecting for underflow, inspired by ANGSD
	if (a > b) {
		return std::log(1 + std::exp(b-a)) + a;
	} // else
	return std::log(1 + std::exp(a-b)) + b;
};

// Can't use Probability type as the sum inside the log will be > 1!
// Also, this function can return a value > 0 which will be normalized at the end
double _protectedSumInLog(double a, double b, double c) {
	// returns log(exp(a)+exp(b)+exp(c)) while protecting for underflow, inspired by ANGSD
	const auto m = std::max({a, b, c});

	return std::log(std::exp(a-m) + std::exp(b-m) + std::exp(c-m)) + m;
};

std::vector<double> chooseLogs(int k) {
	std::vector<double> ret;
	ret.reserve(k + 1);

	// now fill
	for (int j = 0; j <= k; j++) { ret.push_back(coretools::chooseLog(k, j)); }
	return ret;
};

TAlleleCountFile *prepareOutputFile(const std::string & type, std::string filePrefix) {
	// create output file
	TAlleleCountFile *alleleCountFile;
	if (type == "default") {
		filePrefix      = filePrefix + "_alleleCounts.txt.gz";
		alleleCountFile = new TAlleleCountFile(filePrefix);
	} else if (type == "withAlleles") {
		filePrefix      = filePrefix + "_alleleCounts.txt.gz";
		alleleCountFile = new TAlleleCountFileWithAlleles(filePrefix);
	} else if (type == "treemix") {
		filePrefix      = filePrefix + "_treemix_alleleCounts.txt.gz";
		alleleCountFile = new TTreeMixFile(filePrefix);
	} else if (type == "flink") {
		filePrefix      = filePrefix + "_flink_alleleCounts.txt.gz";
		alleleCountFile = new TFlinkFile(filePrefix);
	} else
		UERROR("Unknown output file type '", type, "'! Use 'default', 'withAlleles', 'treemix' or 'flink'.");

	logfile().list("Will write estimated allele counts to file '" + filePrefix + "'.");
	alleleCountFile->openFileToWrite(filePrefix);

	return alleleCountFile;
}

} // namespace impl

//-------------------------------------------------
// TSiteAlleleFrequencyLikelihoods
//-------------------------------------------------
const std::vector<double> &TSiteAlleleFrequencyLikelihoods::_getLogChoose(size_t counts) {
	if (counts >= log_choose.size()) {
		log_choose.resize(counts + 1);
	}
	if (log_choose[counts].empty()) log_choose[counts] = impl::chooseLogs(counts);
	return log_choose[counts];
};

void TSiteAlleleFrequencyLikelihoods::_fillLog(const TSampleLikelihoods * data, uint32_t numSamples) {
	using BG = genometools::BiallelicGenotype;
	// Calculating allele frequency likelihoods according to Nielsen et al. (2012) PLoS One, page 3
	// adapted to also work for haploid individuals (which only have likelihoods for genotypes 0 and 1)
	// all calculations done in log

	// set all h_j = 0
	std::vector<double> alleleFrequencyLikelihoods_h(numSamples*2 + 1, 0.0);
	size_t numAlleleCounts = 0;

	// find first individual  with data
	size_t s = 0;
	while (s < numSamples && data[s].isMissing()) { ++s; }

	if (s < numSamples) {
		// initialize with first individual
		if (data[s].isHaploid()) {
			alleleFrequencyLikelihoods_h[0] = (LogProbability)data[s][BG::haploidFirst];
			alleleFrequencyLikelihoods_h[1] = (LogProbability)data[s][BG::haploidSecond];
			numAlleleCounts                 = 1;
		} else {
			alleleFrequencyLikelihoods_h[0] = (LogProbability)data[s][BG::homoFirst];
			alleleFrequencyLikelihoods_h[1] = (LogProbability)data[s][BG::het] + logOf2;
			alleleFrequencyLikelihoods_h[2] = (LogProbability)data[s][BG::homoSecond];
			numAlleleCounts                 = 2;
		}

		// Recursion
		for (++s; s < numSamples; ++s) {
			if (!data[s].isMissing()) {
				int j = numAlleleCounts;

				if (data[s].isHaploid()) {
					// first fill new ones to avoid multiplication with zero (relevant in log)
					alleleFrequencyLikelihoods_h[j + 1] =
						(LogProbability)data[s][BG::haploidSecond] + alleleFrequencyLikelihoods_h[j];

					// now fill those already used
					for (; j > 0; j--) {
						alleleFrequencyLikelihoods_h[j] = impl::_protectedSumInLog(
							(LogProbability)data[s][BG::haploidSecond] + alleleFrequencyLikelihoods_h[j - 1],
							(LogProbability)data[s][BG::haploidFirst] + alleleFrequencyLikelihoods_h[j]);
					}

					// special case for j=0
					alleleFrequencyLikelihoods_h[0] =
						(LogProbability)data[s][BG::haploidFirst] + alleleFrequencyLikelihoods_h[0];

					// increase total number of haplotypes by one
					numAlleleCounts += 1;
				} else {
					// first fill new ones to avoid multiplication with zero (relevant in log)
					alleleFrequencyLikelihoods_h[j + 2] =
						(LogProbability)data[s][BG::homoSecond] + alleleFrequencyLikelihoods_h[j];
					alleleFrequencyLikelihoods_h[j + 1] = impl::_protectedSumInLog(
						(LogProbability)data[s][BG::homoSecond] + alleleFrequencyLikelihoods_h[j - 1],
						(LogProbability)data[s][BG::het] + alleleFrequencyLikelihoods_h[j] + logOf2);

					// now fill those already used
					for (; j > 1; j--) {
						alleleFrequencyLikelihoods_h[j] = impl::_protectedSumInLog(
							(LogProbability)data[s][BG::homoSecond] + alleleFrequencyLikelihoods_h[j - 2],
							(LogProbability)data[s][BG::het] + alleleFrequencyLikelihoods_h[j - 1] + logOf2,
							(LogProbability)data[s][BG::homoFirst] + alleleFrequencyLikelihoods_h[j]);
					}

					// special case for j=1,0
					alleleFrequencyLikelihoods_h[1] =
						impl::_protectedSumInLog((LogProbability)data[s][BG::het] + alleleFrequencyLikelihoods_h[0] +
											   logOf2, // 2*het does not appear in equation in paper, but must be there
										   (LogProbability)data[s][BG::homoFirst] + alleleFrequencyLikelihoods_h[1]);

					alleleFrequencyLikelihoods_h[0] =
						(LogProbability)data[s][BG::homoFirst] + alleleFrequencyLikelihoods_h[0];

					// increase total number of haplotypes by two
					numAlleleCounts += 2;
				}
			}
		}

		// Termination: add binomial coefficient
		const auto& logChoose = _getLogChoose(numAlleleCounts);
		for (size_t j = 0; j < numAlleleCounts + 1; j++) {
			// numerical accuracy may raely lead to a value very slightly above 0.0. std::min is used to avoid an error
			// when storing as LogProbability.
			log_alleleFrequencyLikelihoods_h.emplace_back(std::min(alleleFrequencyLikelihoods_h[j] - logChoose.at(j), 0.0));
		}
	}
};

void TSiteAlleleFrequencyLikelihoods::_fillNatural(const TSampleLikelihoods * data, uint32_t numSamples) {
	using BG = genometools::BiallelicGenotype;
	// Calculating allele frequency likelihoods according to Nielsen et al. (2012) PLoS One, page 3
	// adapted to also work for haploid individuals (which only have likelihoods for genotypes 0 and 1)
	std::vector<double> alleleFrequencyLikelihoods_h(numSamples*2 + 1, 0.0);

	size_t numAlleleCounts = 0;

	// find first individual  with data
	size_t s = 0;
	while (s < numSamples && data[s].isMissing()) { ++s; }

	if (s < numSamples) {
		// initialize with first individual
		if (data[s].isHaploid()) {
			alleleFrequencyLikelihoods_h[0] = (Probability)data[s][BG::haploidFirst];
			alleleFrequencyLikelihoods_h[1] = (Probability)data[s][BG::haploidSecond];
			numAlleleCounts                 = 1;
		} else {
			alleleFrequencyLikelihoods_h[0] = (Probability)data[s][BG::homoFirst];
			alleleFrequencyLikelihoods_h[1] = (Probability)data[s][BG::het] * 2.0;
			alleleFrequencyLikelihoods_h[2] = (Probability)data[s][BG::homoSecond];
			numAlleleCounts                 = 2;
		}

		// Recursion
		for (++s; s < numSamples; ++s) {
			if (!data[s].isMissing()) {
				int j = numAlleleCounts;

				if (data[s].isHaploid()) {
					// first fill new ones to avoid multiplication with zero (relevant in log, kept here for code
					// consistency)
					alleleFrequencyLikelihoods_h[j + 1] =
						(Probability)data[s][BG::haploidSecond] * alleleFrequencyLikelihoods_h[j];

					// now fill those already used
					for (; j > 0; j--) {
						alleleFrequencyLikelihoods_h[j] =
							(Probability)data[s][BG::haploidSecond] * alleleFrequencyLikelihoods_h[j - 1] +
							(Probability)data[s][BG::haploidFirst] * alleleFrequencyLikelihoods_h[j];
					}

					// special case for j=0
					alleleFrequencyLikelihoods_h[0] =
						(Probability)data[s][BG::haploidFirst] * alleleFrequencyLikelihoods_h[0];

					// increase total number of haplotypes by one
					numAlleleCounts += 1;
				} else {
					// first fill new ones to avoid multiplication with zero (relevant in log, kept here to code
					// consistent)
					alleleFrequencyLikelihoods_h[j + 2] =
						(Probability)data[s][BG::homoSecond] * alleleFrequencyLikelihoods_h[j];
					alleleFrequencyLikelihoods_h[j + 1] =
						(Probability)data[s][BG::homoSecond] * alleleFrequencyLikelihoods_h[j - 1] +
						2.0 * (Probability)data[s][BG::het] * alleleFrequencyLikelihoods_h[j];

					// now fill those already used
					for (; j > 1; j--) {
						alleleFrequencyLikelihoods_h[j] =
							(Probability)data[s][BG::homoSecond] * alleleFrequencyLikelihoods_h[j - 2] +
							2.0 * (Probability)data[s][BG::het] * alleleFrequencyLikelihoods_h[j - 1] +
							(Probability)data[s][BG::homoFirst] * alleleFrequencyLikelihoods_h[j];
					}

					// special case for j=1,0
					alleleFrequencyLikelihoods_h[1] =
						2.0 * (Probability)data[s][BG::het] *
							alleleFrequencyLikelihoods_h[0] // 2*het does not appear in equation in paper, but must be
															// there
						+ (Probability)data[s][BG::homoFirst] * alleleFrequencyLikelihoods_h[1];
					alleleFrequencyLikelihoods_h[0] =
						(Probability)data[s][BG::homoFirst] * alleleFrequencyLikelihoods_h[0];

					// increase total number of haplotypes by two
					numAlleleCounts += 2;
				}
			}
		}

		// Termination: put in log and add binomial coefficient
		const auto &logChoose = _getLogChoose(numAlleleCounts);
		for (size_t j = 0; j < numAlleleCounts + 1; j++) {
			// numerical accuracy may raely lead to a value very slightly above 0.0. std::min is used to avoid an error
			// when storing as LogProbability.
			log_alleleFrequencyLikelihoods_h.emplace_back(std::min(log(alleleFrequencyLikelihoods_h[j]) - logChoose.at(j), 0.0));
		}
	}
};

void TSiteAlleleFrequencyLikelihoods::fill(TSampleLikelihoods *data, uint32_t numSamples, bool resetMissing) {
	// smallest likelihood is 10^-25.5 (phred 255).
	// A double can store up to 10^-308.
	// Hence we can store up to (10^25.5)^12 without underflow
	if (resetMissing) {
		for (size_t i = 0; i < numSamples; ++i) {
			data[i].resetMissing();
		}
	}
	log_alleleFrequencyLikelihoods_h.clear();
	if (numSamples > 12) {
		_fillLog(data, numSamples);
	} else {
		_fillNatural(data, numSamples);
	}
};

void TSiteAlleleFrequencyLikelihoods::write(gz::ogzstream &file) {
	for (const auto& afl: log_alleleFrequencyLikelihoods_h) file << "\t" << afl;
};

size_t TSiteAlleleFrequencyLikelihoods::MLAlleleCount() {
	// return 0 in case of no data
	if (Nalleles() == 0) return 0;

	// first find ML and store all indexes that are at ML
	const auto ML = *std::max_element(log_alleleFrequencyLikelihoods_h.begin(), log_alleleFrequencyLikelihoods_h.end());

	// now store all index at ML
	std::vector<size_t> MLEs;
	for (size_t j = 0; j < log_alleleFrequencyLikelihoods_h.size(); j++) {
		if (log_alleleFrequencyLikelihoods_h[j] == ML) { MLEs.push_back(j); }
	}

	// now choose randomly among those ate MLE
	return MLEs[randomGenerator().sample(MLEs.size())];
};

const std::vector<coretools::LogProbability> &
TSiteAlleleFrequencyLikelihoods::getLogAlleleFrequencyLikelihoods() const {
	return log_alleleFrequencyLikelihoods_h;
};

//-------------------------------------------------
// TAlleleCountEstimator
//-------------------------------------------------
void runCounts() {
	// read samples
	genometools::TPopulationSamples samples;
	if (parameters().exists("samples")) samples.readSamples(parameters().get<std::string>("samples"));

	// open VCF reader
	std::string vcfFilename = parameters().get<std::string>("vcf");
	logfile().startIndent("Reading genotype likelihoods from VCF file '" + vcfFilename + "':");
	genometools::TPopulationLikelihoodReaderLocus reader(false);
	reader.openVCF(vcfFilename);
	logfile().endIndent();

	// Match samples
	if (samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	else
		samples.readSamplesFromVCFNames(reader.getSampleVCFNames());

	// prepare site allele frequency likelihood calculators
	std::vector<TSiteAlleleFrequencyLikelihoods> saf;
	for (size_t p = 0; p < samples.numPopulations(); p++) {
		saf.emplace_back();
	}

	// create out file
	std::string tmp                   = coretools::str::extractBeforeLast(vcfFilename, ".vcf");
	std::string outname               = parameters().get<std::string>("out", tmp);
	std::string type                  = parameters().get<std::string>("outFormat", "default");
	logfile().list("Setting outFormat to " + type + ". (parameter 'outFormat')");
	std::unique_ptr<TAlleleCountFile> alleleCountFile{impl::prepareOutputFile(type, outname)};

	// write header
	alleleCountFile->writeHeader(samples);

	// initialize variables for vcf-file
	struct timeval start;
	gettimeofday(&start, NULL);
	genometools::TPopulationLikehoodLocus<TSampleLikelihoods> data(samples.numSamples());

	// run through VCF file
	logfile().startIndent("Parsing VCF file and estimating allele counts:");
	while (reader.readDataFromVCF(data, samples)) {
		// write chromosome and position
		alleleCountFile->writePosition(reader);

		// print MLE count for each population
		for (size_t p = 0; p < samples.numPopulations(); p++) {
			// calculate allele frequency likelihoods
			saf[p].fill(&data[samples.startIndex(p)], samples.numSamplesInPop(p));

			// and print MLE counts
			alleleCountFile->writeCounts(saf[p].MLAlleleCount(), saf[p].Nalleles(), p);
		}
		alleleCountFile->endl();
	}

	// report final status
	logfile().endIndent();
	reader.concludeFilters();
	logfile().endIndent();
};


void runLikelihoods() {
	// TODO: write proper saf
	// read samples
	genometools::TPopulationSamples samples;
	if (parameters().exists("samples")) samples.readSamples(parameters().get<std::string>("samples"));

	// open VCF reader
	const auto vcfFilename = parameters().get<std::string>("vcf");
	logfile().startIndent("Reading genotype likelihoods from VCF file '" + vcfFilename + "':");
	genometools::TPopulationLikelihoodReaderLocus reader(false);
	reader.openVCF(vcfFilename);
	logfile().endIndent();

	// Match samples
	if (samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	else
		samples.readSamplesFromVCFNames(reader.getSampleVCFNames());

	// prepare site allele frequency likelihood calculators
	const auto tmp     = coretools::str::readBeforeLast(vcfFilename, ".vcf");
	const auto outname = parameters().get("out", tmp);
	logfile().list("Will write estimated allele counts to file '", outname, "'.");

	std::vector<TSiteAlleleFrequencyLikelihoods> saf;
	std::vector<TSafFile> safFiles;
	for (size_t p = 0; p < samples.numPopulations(); p++) {
		saf.emplace_back();
		safFiles.emplace_back(coretools::str::toString(outname, "_", p), samples.numSamplesInPop(p));
	}

	genometools::TPopulationLikehoodLocus<TSampleLikelihoods> data(samples.numSamples());

	// run through VCF file
	logfile().startIndent("Parsing VCF file and estimating allele counts:");
	while (reader.readDataFromVCF(data, samples)) {
		// print MLE count for each population
		for (size_t p = 0; p < samples.numPopulations(); p++) {
			// calculate allele frequency likelihoods
			saf[p].fill(&data[samples.startIndex(p)], samples.numSamplesInPop(p), true);
			safFiles[p].write(reader.chr(), reader.position(), saf[p].getLogAlleleFrequencyLikelihoods());
		}
	}

	// report final status
	logfile().endIndent();
	reader.concludeFilters();
	logfile().endIndent();
};


void runTransform() {
	// initialize variables for vcf-file
	struct timeval start;
	gettimeofday(&start, NULL);

	// get parameters for in and output
	std::string countsFileName = parameters().get<std::string>("transform");
	std::string tmp(coretools::str::readBeforeLast(countsFileName, "_alleleCounts.txt.gz"));
	std::string outname        = parameters().get<std::string>("out", tmp);
	std::string type           = parameters().get<std::string>("outFormat", "default");

	logfile().list("Use option 'outFormat' to produce alleleCounts file in different formats.");

	// open input file	
	TAlleleCountReader file(countsFileName);	

	// create output file
	std::unique_ptr<TAlleleCountFile> alleleCountFile(impl::prepareOutputFile(type, outname));

	// write header
	alleleCountFile->writeHeader(file.populationNames());

	// run through allele count file
	logfile().startIndent("Converting allele counts file:");
	for(const auto& ac : file){		
		// write chromosome and position
		alleleCountFile->writePosition(ac.chr, ac.pos);

		// write counts
		for (unsigned int p = 0; p < ac.numPopulations(); p++) {				
			alleleCountFile->writeCounts(ac[p].minor, ac[p].total, p);
		}
		alleleCountFile->endl();		
	}

	// report final status
	logfile().endIndent();
};

void TAlleleCounter::run() {
	if (parameters().exists("dosaf")) {
		runLikelihoods();
	} else if (parameters().exists("transform")) {
		runTransform();
	} else {
		runCounts();
	}
}

}; // namespace PopulationTools
