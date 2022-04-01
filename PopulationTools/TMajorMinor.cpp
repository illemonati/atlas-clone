/*
 * TMajorMinor.cpp
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#include "TMajorMinor.h"

#include "GenotypeTypes.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "debugtools.h"
#include <algorithm>

namespace PopulationTools {


using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::instances::randomGenerator;
using genometools::AllelicCombination;
using genometools::Base;
using genometools::index;
using GLF::TMultiGLFData;
//---------------------------------------------------
// TMajorMinorEstimatorBase
//---------------------------------------------------

namespace /* anonymous */ {

constexpr auto useAllAlleleicCombinations() {
	std::array<AllelicCombination, index(AllelicCombination::max)> used{};
	for (auto ac = AllelicCombination::min; ac < AllelicCombination::max; ++ac) {
		used[index(ac)] = ac;
	}
	return used;
};

constexpr auto useAllelicCombinationsThatContain(Base base) {
	std::array<AllelicCombination, index(AllelicCombination::max)> used{};

	// get the three alleleic combinations that contain base
	size_t i = 0;
	for (auto b = Base::min; b < Base::max; ++b) {
		if (b != base) { used[i++] = allelicCombination(base, b); }
	}
	return used;
};

} // namespace

void TMajorMinorEstimatorBase::_chooseBestAllelicCombinationAmongThoseWithEqualScores() {
	// identify max L10L
	const auto first_max = std::max_element(_L10L_perCombination.begin(), _L10L_perCombination.end());
	L10L = *first_max;

	// select MLE combination
	std::array<AllelicCombination, index(AllelicCombination::max)> best_combinations;
	size_t i = 0;
	for (auto ac = AllelicCombination(std::distance(_L10L_perCombination.begin(), first_max)); ac < AllelicCombination::max; ++ac) {
		if (_L10L_perCombination[ac] == L10L) best_combinations[i++] = ac;
	}
	bestAllelicCombination = (i == 1 ? best_combinations.front() : best_combinations[randomGenerator().sample(i)]);
};

void TMajorMinorEstimatorBase::_estimateMajorMinor(const TMultiGLFData &data) {
	using coretools::Log10Probability;

	// which one is major?
	if (_genotypeFrequencies.alleleFrequency() < 0.5) {
		major = first(bestAllelicCombination);
		minor = second(bestAllelicCombination);
	} else {
		major = second(bestAllelicCombination);
		minor = first(bestAllelicCombination);

		// also flip frequencies
		_genotypeFrequencies.flip();
	}

	// calculate variant quality
	const auto refHom                  = genometools::genotype(major, major);
	Log10Probability LL_fixed_glfPhred = 0.0;
	for (uint32_t i = 0; i < data.size(); ++i) {
		if (data[i].hasData()) {
			if (data[i].isHaploid())
				LL_fixed_glfPhred += (Log10Probability)data[i][major];
			else
				LL_fixed_glfPhred += (Log10Probability)data[i][refHom];
		}
	}

	// set variant quality
	variantQuality = LL_fixed_glfPhred > L10L ? Log10Probability(0.0) : LL_fixed_glfPhred - L10L;
};

void TMajorMinorEstimatorBase::estimateMajorMinor(const TMultiGLFData &data) {
	constexpr auto used = useAllAlleleicCombinations();
	_findMLAllelicCombination(data, used);
	_estimateMajorMinor(data);
};

void TMajorMinorEstimatorBase::estimateMajorMinor(const TMultiGLFData &data, Base base) {
	const auto used = useAllelicCombinationsThatContain(base);
	_findMLAllelicCombination(data, used);
	_estimateMajorMinor(data);
};

//---------------------------------------------------
// TMajorMinorEstimatorSkotte
//---------------------------------------------------


TMajorMinorEstimatorSkotte::TMajorMinorEstimatorSkotte(double EpsilonF) : epsilonF(EpsilonF){
	using BG = genometools::BiallelicGenotype;

	// diploid
	priorGenotypeFrequencies[BG::homoFirst]  = 0.25;
	priorGenotypeFrequencies[BG::het]        = 0.50;
	priorGenotypeFrequencies[BG::homoSecond] = 0.25;

	// haploid
	priorGenotypeFrequencies[BG::haploidFirst]  = 0.50;
	priorGenotypeFrequencies[BG::haploidSecond] = 0.50;
};

void TMajorMinorEstimatorSkotte::_findMLAllelicCombination(const TMultiGLFData &data, const TAlleleicCombinations& used) {
	// calculate L10L for each allelic combination used
	const auto N = used[4] == AllelicCombination::min ? 4 : used.size(); // Kind of ugly but doing the trick
	for (size_t i = 0; i < N; ++i) {
		const auto ac = used[i];
		fill(_genotypeLikelihoods, data, ac);
		_L10L_perCombination[ac] =
		    priorGenotypeFrequencies.calculateLog10Likelihood(_genotypeLikelihoods, _genotypeLikelihoods.size());
	}

	// pick combination with highest likelihood
	_chooseBestAllelicCombinationAmongThoseWithEqualScores();

	// now estimate genotype frequencies at MLE allelic combination
	fill(_genotypeLikelihoods, data, bestAllelicCombination);
	_genotypeFrequencies.estimate(_genotypeLikelihoods, _genotypeLikelihoods.size(), epsilonF);

	// calculate likelihood again with better genotype frequencies
	_L10L_perCombination[bestAllelicCombination] =
	    _genotypeFrequencies.calculateLog10Likelihood(_genotypeLikelihoods, _genotypeLikelihoods.size());
	L10L = _L10L_perCombination[bestAllelicCombination];
};

//---------------------------------------------------
// TMajorMinorEstimatorMLE
//---------------------------------------------------
coretools::Log10Probability TMajorMinorEstimatorMLE::estimateGenotypeFrequencies(const TMultiGLFData &data,
										 AllelicCombination ac) {
	using genometools::index;
	fill(_genotypeLikelihoods, data, ac);
	const auto sz = _genotypeLikelihoods.size();
	tmpGenotypeFrequencies[index(ac)].estimate(_genotypeLikelihoods, sz, epsilonF);
	return tmpGenotypeFrequencies[index(ac)].calculateLog10Likelihood(_genotypeLikelihoods, sz);
};

void TMajorMinorEstimatorMLE::_findMLAllelicCombination(const TMultiGLFData &data, const TAlleleicCombinations& used) {
	// calculate L10L for each allelic combination
	const auto N = used[4] == AllelicCombination::min ? 4 : used.size(); // Kind of ugly but doing the trick
	for (size_t i = 0; i < N; ++i) {
		const auto ac = used[i];
		_L10L_perCombination[ac] = estimateGenotypeFrequencies(data, ac);
	}

	// pick combination
	_chooseBestAllelicCombinationAmongThoseWithEqualScores();
	_genotypeFrequencies.set(tmpGenotypeFrequencies[index(bestAllelicCombination)]);
};

//---------------------------------------------------
// TMajorMinor
//---------------------------------------------------
void estimateMajorMinor() {
	// open GLF files
	GLF::TGlfMultiReader glfReader;
	glfReader.openGLFs();
	glfReader.setAllActive();

	bool hasReference = false;

	// add reference, if provided
	if (parameters().parameterExists("fasta")) {
		logfile().list("Will only identify the most likely alternative allele (argument: fasta)");
		std::string fastaFile = parameters().getParameter<std::string>("fasta");
		logfile().list("Reading reference sequence from '" + fastaFile + "'");
		glfReader.addReference(fastaFile);
		hasReference = true;
	}

	// estimation method
	std::string method = parameters().getParameterWithDefault<std::string>("method", "MLE");
	TMajorMinorEstimatorBase *MMEstimator;
	double maxF = parameters().getParameterWithDefault("maxF", 0.0000001);
	if (method == "Skotte") {
		logfile().list("Will estimate major / minor alleles using the Skotte method with maxF ", maxF,
			      ". (parameters method and maxF)");
		MMEstimator = new TMajorMinorEstimatorSkotte(maxF);
	} else if (method == "MLE") {
		logfile().list("Will estimate major / minor alleles using the MLE method with maxF ", maxF,
			      ". (parameters method and maxF)");
		MMEstimator = new TMajorMinorEstimatorMLE(maxF);
	} else
		throw "Unknown MajorMinor method '" + method + "'!";

	bool usePhredLikelihoods = parameters().parameterExists("phredLik");
	if (usePhredLikelihoods) {
		logfile().list("Will write phred-scaled likelihoods. (parameter phredLik)");
	} else {
		logfile().list("Will write raw likelihoods. (use phredLik to phred-scale)");
	}

	// read filters
	uint32_t minSamplesWithData = 1;
	genometools::PhredIntProbability minVariantQuality{0};
	if (parameters().parameterExists("printAll")) {
		logfile().list("Will all sites and samples. (parameter printAll)");
		minSamplesWithData = 0;
		minVariantQuality  = 0;
	} else {
		minSamplesWithData = parameters().getParameterWithDefault<uint32_t>("minSamplesWithData", 1);
		if (minSamplesWithData > 0) {
			logfile().list("Will only print sites for which at least ", minSamplesWithData,
				      " samples have data. (parameter minSamplesWithData)");
		}

		minVariantQuality = parameters().getParameterWithDefault<genometools::PhredIntProbability>(
		    "minVariantQual", genometools::PhredIntProbability::highest());
		if (minVariantQuality > genometools::PhredIntProbability::highest()) {
			logfile().list("Will only print sites with variant quality >= ", minVariantQuality,
				      " samples have data. (parameter minVariantQual)");
		}
	}

	if (minSamplesWithData > 0) {
		glfReader.onlyJumpToPositionsWithData(true);
	} else {
		glfReader.onlyJumpToPositionsWithData(false);
	}

	// limit input
	long limitSites = parameters().getParameterWithDefault("limitSites", 0);
	if (limitSites > 0) logfile().list("Will stop at input position ", limitSites, ". (parameter 'limitSites')");
	if (limitSites < 0) throw "maxPos cannot be negative!";

	// filename tag
	std::string outname = parameters().getParameterWithDefault<std::string>("out", "ATLAS_majorMinor");
	logfile().list("Will write output files with tag '" + outname + "'. (parameter 'out')");

	// get sample names. Make sure they are unique in the vcf
	std::vector<std::string> sampleNames;
	if (parameters().parameterExists("sampleNames")) {
		logfile().startIndent("Using the following sample names (parameter 'sampleNames'):");
		parameters().fillParameterIntoContainer("sampleNames", sampleNames, ',');

		if (sampleNames.size() != glfReader.numActiveSamples()) {
			throw "Number of provided same names does not match number of GLF files!";
		}

		// report
		auto glfNames = glfReader.namesOfActiveFiles();
		for (uint32_t i = 0; i < glfReader.numActiveSamples(); ++i) {
			logfile().list(glfNames[i], " -> ", sampleNames[i]);
		}
		logfile().endIndent();
	} else {
		logfile().list("Will deduce sample names from GLF file names. (use 'sampleNames' to provide alternative names)");
		sampleNames = glfReader.sampleNamesOfActiveFiles();

		// if there are duplicates, add suffix
		bool foundDuplicates = false;
		for (auto &s : sampleNames) {
			uint32_t c = 0;
			for (auto &t : sampleNames) {
				if (t == s) {
					++c;
					if (c > 1) {
						t = t + "." + coretools::str::toString(c);

						// report
						if (!foundDuplicates) {
							logfile().startIndent("Duplicate samples will be rename as follows:");
							foundDuplicates = true;
						}
						logfile().list(s, " -> ", t);
					}
				}
			}
		}
		if (foundDuplicates) { logfile().endIndent(); }
	}

	// open vcf file
	GLF::TGlfMultiReaderVcf vcf(outname + ".vcf.gz", "ATLAS_GLF_Caller", sampleNames, usePhredLikelihoods);

	// vars
	logfile().startIndent("Parsing through glf files:");
	coretools::TTimer timer;
	long counter = 0;

	while (glfReader.readNext()) {
		++counter;
		// filter on missingness
		if (glfReader.numActiveSamplesWithData() >= minSamplesWithData) {
			// do major minor
			if (hasReference) {
				Base ref = glfReader.refBase();
				MMEstimator->estimateMajorMinor(glfReader.data, ref);

				// write to VCF
				if (MMEstimator->variantQuality >= minVariantQuality) {
					if (MMEstimator->major == ref) {
						vcf.writeSite(glfReader.chr(), glfReader.position(), MMEstimator->variantQuality,
							      glfReader.data, MMEstimator->major, MMEstimator->minor);
					} else {
						vcf.writeSite(glfReader.chr(), glfReader.position(), MMEstimator->variantQuality,
							      glfReader.data, MMEstimator->minor, MMEstimator->major);
					}
				}
			} else {
				MMEstimator->estimateMajorMinor(glfReader.data);

				// write to VCF
				if (MMEstimator->variantQuality >= minVariantQuality) {
					vcf.writeSite(glfReader.chr(), glfReader.position(), MMEstimator->variantQuality, glfReader.data,
					              MMEstimator->major, MMEstimator->minor);
				}
			}
		} // end filter on missingness

		// report progress
		if (counter % 1000000 == 0) {
			logfile().list("Parsed ", counter, " positions in ", timer.formattedTime(), " min.");
		}

		// break?
		if (limitSites > 0 && counter == limitSites) break;
	}

	logfile().list("Reached end of glf files!");
	logfile().removeIndent();

	// clean storage
	delete MMEstimator;
};

}; // namespace PopulationTools
