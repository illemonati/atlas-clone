/*
 * TMajorMinor.cpp
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#include "TMajorMinor.h"

#include <algorithm>
#include <array>
#include <exception>
#include <memory>
#include <set>
#include <stddef.h>
#include <stdint.h>
#include <string>
#include <tuple>
#include <vector>

#include "TGlfMultiReader.h"
#include "coretools/Containers/TDualArray.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Containers/TView.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Main/TTask.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/TTimer.h"
#include "coretools/Types/probability.h"
#include "coretools/Types/strongTypes.h"
#include "coretools/Types/weakTypes.h"
#include "coretools/devtools.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "genometools/TGenotypeFrequencies.h"

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

namespace impl {

using TAlleleicCombinationData = coretools::TStrongArray<coretools::Log10Probability, genometools::AllelicCombination>;

constexpr coretools::TDualArray<AllelicCombination, 3, index(AllelicCombination::max)>
useAllelicCombinationsThatContain(Base base) {
	using AC = AllelicCombination;
	switch (base) {
	case Base::A: return std::array{AC::AC, AC::AG, AC::AT};
	case Base::C: return std::array{AC::AC, AC::CG, AC::CT};
	case Base::G: return std::array{AC::AG, AC::CG, AC::GT};
	case Base::T: return std::array{AC::AT, AC::CT, AC::GT};
	case Base::N: return std::array{AC::AC, AC::AG, AC::AT, AC::CG, AC::CT, AC::GT};
	}
};

AllelicCombination chooseBestAllelicCombination(const TAlleleicCombinationData &acd) {
	return randomGenerator()
		.sampleIndexOfMaxima<TAlleleicCombinationData, AllelicCombination, index(AllelicCombination::max)>(acd);
};

} // namespace impl

struct TMMData {
	coretools::Probability MAF;
	genometools::PhredIntProbability variantQuality;
	genometools::Base major;
	genometools::Base minor;
};

template<typename AllelicFinder>
static TMMData estimate(coretools::TConstView<GLF::TMultiGLFDataSample> data, double maxF, genometools::Base base = genometools::Base::N) {
	using coretools::Log10Probability;

	auto [_bestAllelicCombination, L, genotypeFrequencies] = AllelicFinder::find(data, maxF, base);
	Base _major, _minor;

	// which one is major?
	if (genotypeFrequencies.alleleFrequency() < 0.5) {
		_major = first(_bestAllelicCombination);
		_minor = second(_bestAllelicCombination);
	} else {
		_major = second(_bestAllelicCombination);
		_minor = first(_bestAllelicCombination);

		// also flip frequencies
		genotypeFrequencies.flip();
	}

	if (_minor == base) { // cannot happen if base == ref
		_minor = _major;
		_major = base;
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
	genometools::PhredIntProbability _variantQuality{LL_fixed_glfPhred > L ? Log10Probability(0.0)
																		   : Log10Probability(LL_fixed_glfPhred - L)};
	return {genotypeFrequencies.MAF(), _variantQuality, _major, _minor};
}

struct TSkotte {
	static auto find(coretools::TConstView<GLF::TMultiGLFDataSample> data, double maxF, genometools::Base base) {
		// calculate L10L for each allelic combination used
		const auto used = impl::useAllelicCombinationsThatContain(base);
		impl::TAlleleicCombinationData Ls{};
		for (auto ac : used) Ls[ac] = 0.;
		for (const auto &d : data) {
			if (!d.hasData()) continue;
			if (d.isHaploid()) {
				for (auto ac : used) {
					Ls[ac] += log10(0.5 * (coretools::Probability)d[genometools::first(ac)] +
									0.5 * (coretools::Probability)d[genometools::second(ac)]);
				}
			} else {
				for (auto ac : used) {
					Ls[ac] += log10(0.25 * (coretools::Probability)d[genometools::homoFirst(ac)] +
									0.5 * (coretools::Probability)d[genometools::het(ac)] +
									0.25 * (coretools::Probability)d[genometools::homoSecond(ac)]);
				}
			}
		}
		AllelicCombination bestAC = impl::chooseBestAllelicCombination(Ls);

		// now estimate genotype frequencies at MLE allelic combination
		auto GLs = fill(data, bestAC);
		genometools::TGenotypeFrequencies GFs;
		GFs.estimate<false>(GLs, GLs.size(), maxF);

		return std::make_tuple(bestAC, GFs.calculateLog10Likelihood(GLs, GLs.size()), GFs);
	}
};

struct TMLE {
	static auto find(coretools::TConstView<GLF::TMultiGLFDataSample> data, double maxF, genometools::Base base) {
		// calculate L10L for each allelic combination
		const auto used = impl::useAllelicCombinationsThatContain(base);

		genometools::TGenotypeFrequencies bestFreqs;
		coretools::Log10Probability bestL = coretools::Log10Probability::lowest();
		AllelicCombination bestAC         = AllelicCombination::min;

		for (const auto ac : used) {
			auto GLs = fill(data, ac);
			genometools::TGenotypeFrequencies freqs;
			freqs.estimate<false>(GLs, GLs.size(), maxF);
			auto L = freqs.calculateLog10Likelihood(GLs, GLs.size());
			if ((L > bestL) || (L == bestL && randomGenerator().getRand() > 0.5)) {
				bestL     = L;
				bestFreqs = freqs;
				bestAC    = ac;
			}
		}
		return std::make_tuple(bestAC, bestL, bestFreqs);
	}
};

struct TApproxL {
	static auto find(coretools::TConstView<GLF::TMultiGLFDataSample> data, double maxF, genometools::Base base) {
		// calculate L10L for each allelic combination
		const auto used = impl::useAllelicCombinationsThatContain(base);

		coretools::Log10Probability bestL = coretools::Log10Probability::lowest();
		AllelicCombination bestAC         = AllelicCombination::min;
		for (const auto ac : used) {
			double L = 0.;
			for (const auto &d: data) {
				if (!d.hasData()) continue;
				if (d.isHaploid()) {
					L+= d[genometools::first(ac)].get() + d[genometools::second(ac)].get();
				} else {
					L+= d[genometools::homoFirst(ac)].get() + d[genometools::het(ac)].get() + d[genometools::homoSecond(ac)].get();
				}
			}
			if ((L > bestL) || (L == bestL && randomGenerator().getRand() > 0.5)) {
				bestL     = L;
				bestAC    = ac;
			}
		}

		auto GLs = fill(data, bestAC);
		genometools::TGenotypeFrequencies GFs;
		GFs.estimate<false>(GLs, GLs.size(), maxF);

		return std::make_tuple(bestAC, GFs.calculateLog10Likelihood(GLs, GLs.size()), GFs);
	}
};

template<typename Estimator> void iterate(double maxF) {
	// open GLF files
	GLF::TGlfMultiReader glfReader;
	glfReader.openGLFs();
	glfReader.setAllActive();

	// add reference, if provided
	if (parameters().parameterExists("fasta")) {
		logfile().list("Will only identify the most likely alternative allele (argument: fasta)");
		const std::string fastaFile = parameters().getParameter<std::string>("fasta");
		logfile().list("Reading reference sequence from '" + fastaFile + "'");
		glfReader.addReference(fastaFile);
	}

	// estimation method
	const std::string method = parameters().getParameterWithDefault<std::string>("method", "MLE");

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
	size_t nMAFMAF                = 0;
	if (minMAF > 0.0) {
		logfile().list("Will filter on a minor allele frequency of ", minMAF, ". (parameter 'minMAF')");
	} else {
		logfile().list("Will keep sites regardless of their minor allele frequency. (use 'minMAF' to filter)");
	}

	// limit input
	const size_t limitSites = parameters().getParameterWithDefault("limitSites", 0);
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
		logfile().list(
			"Will deduce sample names from GLF file names. (use 'sampleNames' to provide alternative names)");
		sampleNames = glfReader.sampleNamesOfActiveFiles();

		// if there are duplicates, add suffix
		bool foundDuplicates = false;
		for (size_t i = 0; i < sampleNames.size(); ++i) {
			for (size_t j = i + 1; j < sampleNames.size(); ++j) {
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
	size_t maxThreads = coretools::instances::parameters().getParameterWithDefault("maxThreads", omp_get_max_threads());
	coretools::instances::logfile().list("Running in parallel with a maximum of ", maxThreads,
										 " threads (argument 'maxThreads')");
#else
	coretools::instances::logfile().list("Not running in parallel");
#endif

	// open vcf file
	GLF::TGlfMultiReaderVcf vcf(outname + ".vcf.gz", "ATLAS_GLF_Caller", sampleNames, usePhredLikelihoods);

	// vars
	logfile().startIndent("Parsing through glf files:");
	coretools::TTimer timer;
	constexpr size_t dCounter = 1000000;
	size_t counter            = 0;
	size_t nextPrint          = dCounter;

	for (auto ids = glfReader.readWindow(); !ids.empty(); ids = glfReader.readWindow()) {
		std::vector<TMMData> data(ids.back() + 1);
#pragma omp parallel for num_threads(maxThreads)
		for (auto i : ids) {
			const Base ref = glfReader.refBase(i); // can be N
			data[i]        = estimate<Estimator>(glfReader.data(i), maxF, ref);
		}

		// pass filter?
		for (auto i : ids) {
			const auto &di = data[i];
			if (di.MAF < minMAF) {
				++nMAFMAF;
				continue;
			}
			if (di.variantQuality < minVariantQuality) {
				++nVariantQuality;
				continue;
			}

			// write to VCF
			vcf.writeSite(glfReader.chr(), glfReader.position(i), di.variantQuality, glfReader.data(i), di.major,
						  di.minor);
		}

		counter += ids.back() + 1;

		// report progress
		if (counter >= nextPrint) {
			logfile().list("Parsed ", nextPrint, " positions in ", timer.formattedTime(), ".");
			while (nextPrint <= counter) nextPrint += dCounter;
		}

		if (limitSites > 0 && counter == limitSites) break;
	}

	logfile().list("Reached end of glf files!");

	if (minVariantQuality > genometools::PhredIntProbability::highest()) {
		logfile().conclude("Filtered ", nVariantQuality, " sites with variant quality < ", minVariantQuality, ".");
	}
	if (minMAF > 0) { logfile().conclude("Filtered ", nMAFMAF, " sites with MAF < ", minMAF, "."); }
	logfile().removeIndent();
};

//---------------------------------------------------
// TMajorMinor
//---------------------------------------------------
void TMajorMinor::run() {
	const std::string method = parameters().getParameterWithDefault<std::string>("method", "MLE");

	const double maxF = parameters().getParameterWithDefault("maxF", 0.0000001);
	if (method == "Skotte") {
		logfile().list("Will estimate major / minor alleles using the Skotte method with maxF ", maxF,
					   ". (parameters method and maxF)");
		iterate<TSkotte>(maxF);
	} else if (method == "MLE") {
		logfile().list("Will estimate major / minor alleles using the MLE method with maxF ", maxF,
					   ". (parameters method and maxF)");
		iterate<TMLE>(maxF);
	} else if (method == "ApproxL") {
		logfile().list("Will estimate major / minor alleles using the MLE method with maxF ", maxF,
					   ". (parameters method and maxF)");
		iterate<TApproxL>(maxF);
	} else {
		UERROR("Unknown MajorMinor method '", method, "'!");
	}
}

}; // namespace PopulationTools
