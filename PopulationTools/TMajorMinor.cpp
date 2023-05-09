/*
 * TMajorMinor.cpp
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#include "TMajorMinor.h"

#include <stddef.h>
#include <stdint.h>
#include <algorithm>
#include <exception>
#include <memory>
#include <vector>

#include "coretools/devtools.h"
#include "genometools/GenotypeTypes.h"
#include "coretools/Containers/TDualArray.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/TTimer.h"
#include "coretools/Types/probability.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Types/strongTypes.h"
#include "coretools/Types/weakTypes.h"


#ifdef _OPENMP
#include "omp.h"
#endif


namespace PopulationTools {


using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::instances::randomGenerator;
using genometools::AllelicCombination;
using genometools::Base;
using coretools::index;
using GLF::TMultiGLFData;
//---------------------------------------------------
// TMajorMinorEstimatorBase
//---------------------------------------------------

namespace /* anonymous */ {
constexpr coretools::TDualArray<AllelicCombination, 3, index(AllelicCombination::max)> useAllelicCombinationsThatContain(Base base) {
	using AC = AllelicCombination;
	switch (base) {
	case Base::A: return std::array{AC::AC, AC::AG, AC::AT};
	case Base::C: return std::array{AC::AC, AC::CG, AC::CT};
	case Base::G: return std::array{AC::AG, AC::CG, AC::GT};
	case Base::T: return std::array{AC::AT, AC::CT, AC::GT};
	case Base::N: return std::array{AC::AC, AC::AG, AC::AT, AC::CG, AC::CT,AC::GT};
	}
};

AllelicCombination chooseBestAllelicCombination(const TAlleleicCombinationData& acd) {
	return randomGenerator().sampleIndexOfMaxima<TAlleleicCombinationData, AllelicCombination, index(AllelicCombination::max)>(acd);
};

} // namespace


void TMajorMinorEstimatorBase::estimateMajorMinor(const TMultiGLFData &data, Base base) {
	using coretools::Log10Probability;

	_findMLAllelicCombination(data, base);

	// which one is major?
	if (_genotypeFrequencies.alleleFrequency() < 0.5) {
		_major = first(_bestAllelicCombination);
		_minor = second(_bestAllelicCombination);
	} else {
		_major = second(_bestAllelicCombination);
		_minor = first(_bestAllelicCombination);

		// also flip frequencies
		_genotypeFrequencies.flip();
	}

	// calculate variant quality
	const auto refHom                  = genometools::genotype(_major, _major);
	Log10Probability LL_fixed_glfPhred = 0.0;
	for (size_t i = 0; i < data.size(); ++i) {
		if (data[i].hasData()) {
			if (data[i].isHaploid())
				LL_fixed_glfPhred += (Log10Probability)data[i][_major];
			else
				LL_fixed_glfPhred += (Log10Probability)data[i][refHom];
		}
	}
	_variantQuality = LL_fixed_glfPhred > _L10L_perCombination[_bestAllelicCombination] ? Log10Probability(0.0) : Log10Probability(LL_fixed_glfPhred - _L10L_perCombination[_bestAllelicCombination]);
};

//---------------------------------------------------
// TMajorMinorEstimatorSkotte
//---------------------------------------------------


TMajorMinorEstimatorSkotte::TMajorMinorEstimatorSkotte(double EpsilonF) : _epsilonF(EpsilonF){
	using BG = genometools::BiallelicGenotype;

	// diploid
	_priorGenotypeFrequencies[BG::homoFirst]  = 0.25;
	_priorGenotypeFrequencies[BG::het]        = 0.50;
	_priorGenotypeFrequencies[BG::homoSecond] = 0.25;

	// haploid
	_priorGenotypeFrequencies[BG::haploidFirst]  = 0.50;
	_priorGenotypeFrequencies[BG::haploidSecond] = 0.50;
};

void TMajorMinorEstimatorSkotte::_findMLAllelicCombination(const TMultiGLFData &data, Base base) {
	// calculate L10L for each allelic combination used
	const auto used = useAllelicCombinationsThatContain(base);
	for (const auto ac: used) {
		fill(_genotypeLikelihoods, data, ac);
		_L10L_perCombination[ac] =
		    _priorGenotypeFrequencies.calculateLog10Likelihood(_genotypeLikelihoods, _genotypeLikelihoods.size());
	}

	// pick combination with highest likelihood
	_bestAllelicCombination = chooseBestAllelicCombination(_L10L_perCombination);

	// now estimate genotype frequencies at MLE allelic combination
	fill(_genotypeLikelihoods, data, _bestAllelicCombination);
	_genotypeFrequencies.estimate<false>(_genotypeLikelihoods, _genotypeLikelihoods.size(), _epsilonF);

	// calculate likelihood again with better genotype frequencies
	_L10L_perCombination[_bestAllelicCombination] =
	    _genotypeFrequencies.calculateLog10Likelihood(_genotypeLikelihoods, _genotypeLikelihoods.size());
};

//---------------------------------------------------
// TMajorMinorEstimatorMLE
//---------------------------------------------------
coretools::Log10Probability TMajorMinorEstimatorMLE::_estimateGenotypeFrequencies(const TMultiGLFData &data,
										 AllelicCombination ac) {
	using coretools::index;
	fill(_genotypeLikelihoods, data, ac);
	const auto sz = _genotypeLikelihoods.size();
	_tmpGenotypeFrequencies[index(ac)].estimate<false>(_genotypeLikelihoods, sz, _epsilonF);
	return _tmpGenotypeFrequencies[index(ac)].calculateLog10Likelihood(_genotypeLikelihoods, sz);
};

void TMajorMinorEstimatorMLE::_findMLAllelicCombination(const TMultiGLFData &data, Base base) {
	// calculate L10L for each allelic combination
	const auto used = useAllelicCombinationsThatContain(base);
	for (const auto ac: used) {
		_L10L_perCombination[ac] = _estimateGenotypeFrequencies(data, ac);
	}

	// pick combination
	_bestAllelicCombination = chooseBestAllelicCombination(_L10L_perCombination);
	_genotypeFrequencies.set(_tmpGenotypeFrequencies[index(_bestAllelicCombination)]);
};

//---------------------------------------------------
// TMajorMinor
//---------------------------------------------------
void TMajorMinor::run() {
	// open GLF files
	GLF::TGlfMultiReader glfReader;
	glfReader.openGLFs();
	glfReader.setAllActive();

	bool hasReference = false;

	// add reference, if provided
	if (parameters().parameterExists("fasta")) {
		logfile().list("Will only identify the most likely alternative allele (argument: fasta)");
		const std::string fastaFile = parameters().getParameter<std::string>("fasta");
		logfile().list("Reading reference sequence from '" + fastaFile + "'");
		glfReader.addReference(fastaFile);
		hasReference = true;
	}

	// estimation method
	const std::string method = parameters().getParameterWithDefault<std::string>("method", "MLE");
	const size_t windowSize  = parameters().getParameterWithDefault<size_t>("window", 64);
	std::vector<std::unique_ptr<TMajorMinorEstimatorBase>> MMEstimator;
	const double maxF = parameters().getParameterWithDefault("maxF", 0.0000001);
	if (method == "Skotte") {
		logfile().list("Will estimate major / minor alleles using the Skotte method with maxF ", maxF,
			      ". (parameters method and maxF)");
		for (size_t i = 0; i < windowSize; ++i) {
			MMEstimator.push_back(std::make_unique<TMajorMinorEstimatorSkotte>(maxF));
		}
	} else if (method == "MLE") {
		logfile().list("Will estimate major / minor alleles using the MLE method with maxF ", maxF,
			      ". (parameters method and maxF)");
		for (size_t i = 0; i < windowSize; ++i) {
			MMEstimator.push_back(std::make_unique<TMajorMinorEstimatorMLE>(maxF));
		}
	} else
		UERROR("Unknown MajorMinor method '", method, "'!");

	const bool usePhredLikelihoods = parameters().parameterExists("phredLik");
	if (usePhredLikelihoods) {
		logfile().list("Will write phred-scaled likelihoods. (parameter phredLik)");
	} else {
		logfile().list("Will write raw likelihoods. (use phredLik to phred-scale)");
	}

	// read filters
	size_t minSamplesWithData = 1;
	genometools::PhredIntProbability minVariantQuality{0};
	size_t nVariantQuality = 0;
	if (parameters().parameterExists("printAll")) {
		logfile().list("Will all sites and samples. (parameter printAll)");
		minSamplesWithData = 0;
		minVariantQuality  = 0;
	} else {
		minSamplesWithData = parameters().getParameterWithDefault<size_t>("minSamplesWithData", 1);
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
	glfReader.minSamplesWithData(minSamplesWithData);

	coretools::Probability minMAF = parameters().getParameterWithDefault("minMAF", 0.0);
	size_t nMAFMAF = 0;
	if (minMAF > 0.0) {
		logfile().list("Will filter on a minor allele frequency of ", minMAF, ". (parameter 'minMAF')");
	} else {
		logfile().list("Will keep sites regardless of their minor allele frequency. (use 'minMAF' to filter)");
	}

	// limit input
	const long limitSites = parameters().getParameterWithDefault("limitSites", 0);
	if (limitSites > 0) logfile().list("Will stop at input position ", limitSites, ". (parameter 'limitSites')");
	if (limitSites < 0) UERROR("maxPos cannot be negative!");

	// filename tag
	const std::string outname = parameters().getParameterWithDefault<std::string>("out", "ATLAS_majorMinor");
	logfile().list("Will write output files with tag '" + outname + "'. (parameter 'out')");

	// get sample names. Make sure they are unique in the vcf
	std::vector<std::string> sampleNames;
	if (parameters().parameterExists("sampleNames")) {
		logfile().startIndent("Using the following sample names (parameter 'sampleNames'):");
		parameters().fillParameterIntoContainer("sampleNames", sampleNames, ',');

		if (sampleNames.size() != glfReader.numActiveSamples()) {
			UERROR("Number of provided same names does not match number of GLF files!");
		}

		// report
		auto glfNames = glfReader.namesOfActiveFiles();
		for (size_t i = 0; i < glfReader.numActiveSamples(); ++i) {
			logfile().list(glfNames[i], " -> ", sampleNames[i]);
		}
		logfile().endIndent();
	} else {
		logfile().list("Will deduce sample names from GLF file names. (use 'sampleNames' to provide alternative names)");
		sampleNames = glfReader.sampleNamesOfActiveFiles();

		// if there are duplicates, add suffix
		bool foundDuplicates = false;
		for (size_t i = 0; i < sampleNames.size(); ++i) {
			for (size_t j = i+1; j < sampleNames.size(); ++j) {
				int counter = 1;
				if (sampleNames[i] == sampleNames[j]) {
					sampleNames[j] += "." + coretools::str::toString(counter++);
					if (!foundDuplicates) {
						logfile().startIndent("Duplicate samples will be rename as follows:");
						foundDuplicates = true;
					}
					logfile().list(sampleNames[i], " -> ", sampleNames[j]);
				}
			}
		}
		if (foundDuplicates) { logfile().endIndent(); }
	}

#ifdef _OPENMP
	size_t maxNumThreads =
		coretools::instances::parameters().getParameterWithDefault("maxNumThreads", omp_get_max_threads());
	coretools::instances::logfile().list("Running in parallel with a maximum of ", maxNumThreads,
										 " threads (argument 'maxNumThreads')");
#else
	coretools::instances::logfile().list("Not running in parallel");
#endif

	// open vcf file
	GLF::TGlfMultiReaderVcf vcf(outname + ".vcf.gz", "ATLAS_GLF_Caller", sampleNames, usePhredLikelihoods);

	// vars
	logfile().startIndent("Parsing through glf files:");
	coretools::TTimer timer;
	long counter = 0;


	for (size_t N = glfReader.readWindow(); N > 0; N = glfReader.readWindow()) {
#pragma omp parallel for num_threads(maxNumThreads)
		for (size_t i = 0; i < N; ++i) {
			if (glfReader.numActiveSamplesWithData(i) < minSamplesWithData) continue;
			const Base ref = glfReader.refBase(i); // can be N
			MMEstimator[i]->estimateMajorMinor(glfReader.data(i), ref);
		}

		// pass filter?
		for (size_t i = 0; i < N; ++i) {
			if (glfReader.numActiveSamplesWithData(i) < minSamplesWithData) continue;
			const Base ref = glfReader.refBase(i); // can be N
			if (MMEstimator[i]->genotypeFrequencies().MAF() < minMAF) {
				++nMAFMAF;
				continue;
			}
			if (MMEstimator[i]->variantQuality() < minVariantQuality) {
				++nVariantQuality;
				continue;
			}

			// write to VCF
			if (hasReference && MMEstimator[i]->minor() == ref) {
				vcf.writeSite(glfReader.chr(), glfReader.position(i), MMEstimator[i]->variantQuality(), glfReader.data(i),
							  ref, MMEstimator[i]->major());
			} else {
				vcf.writeSite(glfReader.chr(), glfReader.position(i), MMEstimator[i]->variantQuality(), glfReader.data(i),
							  MMEstimator[i]->major(), MMEstimator[i]->minor());
			}
		}

		counter += N;

		// report progress
		if (counter % 100000 == 0) {
			logfile().list("Parsed ", counter, " positions in ", timer.formattedTime(), ".");
		}

		if (limitSites > 0 && counter == limitSites) break;
	}

	logfile().list("Reached end of glf files!");

	if (minVariantQuality > genometools::PhredIntProbability::highest()) {
		logfile().conclude("Filtered ", nVariantQuality, " sites with variant quality < ", minVariantQuality, ".");
	}
	if (minMAF > 0) {
		logfile().conclude("Filtered ", nMAFMAF, " sites with MAF < ", minMAF, ".");
	}
	logfile().removeIndent();
};

}; // namespace PopulationTools
