/*
 * TMajorMinor.cpp
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#include "TMajorMinor.h"
#include "TMajorMinorMLE.h"
#include "TSkotte.h"

#include <algorithm>

#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Types/probability.h"

#include "genometools/GLF/GLF.h"
#include "genometools/GLF/TMultiGLFTraverser.h"
#include "genometools/Genotypes/TFrequencies.h"
#include "genometools/TAlleles.h"
#include "genometools/TFastaReader.h"
#include "genometools/VCF/TVCFWriter.h"

#include "TBgzWriter.h"

namespace PopulationTools {

using coretools::Log10Probability;
using coretools::P;
using coretools::Probability;
using coretools::PhredInt;
using coretools::TConstView;
using coretools::instances::logfile;
using coretools::instances::parameters;
using genometools::Base;
using genometools::TGLFEntry;

namespace MajorMinor {
coretools::Log10Probability LLFixedAllele(coretools::TConstView<genometools::TGLFEntry> data, genometools::Base major) {
	using coretools::Log10Probability;
	const auto refHom = genometools::genotype(major, major);
	Log10Probability LL_fixed{0.0};
	for (size_t i = 0; i < data.size(); ++i) {
		if (data[i].depth) {
			if (data[i].isHaploid())
				LL_fixed += (Log10Probability)data[i][major];
			else
				LL_fixed += (Log10Probability)data[i][refHom];
		}
	}
	return LL_fixed;
}

std::array<genometools::AllelicCombination, 3> useAllelicCombinationsThatContain(genometools::Base base) {
	using genometools::Base;
	DEBUG_ASSERT(base != Base::N);
	using AC = genometools::AllelicCombination;
	switch (base) {
	case Base::A: return std::array{AC::AC, AC::AG, AC::AT};
	case Base::C: return std::array{AC::AC, AC::CG, AC::CT};
	case Base::G: return std::array{AC::AG, AC::CG, AC::GT};
	default: return std::array{AC::AT, AC::CT, AC::GT};
	}
};
} // namespace MajorMinor

template<typename Estimator> void iterate(double maxF) {
	// open GLF files
	genometools::TMultiGLFTraverser multiTraverser;
	genometools::TFastaReader fastaReader;

	// use known alleles or reference allele, if provided

	bool filterN = false;

	if (multiTraverser.alleles()) {
		logfile().list("Will limit analysis to sites with known alleles (parameter 'alleles'):");
	} else if (parameters().exists("fasta")) {
		logfile().list(
		    "Will use reference allele and only identify the most likely alternative allele. (argument: fasta)");
		const std::string fastaFile = parameters().get<std::string>("fasta");
		fastaReader.open(fastaFile);
		filterN = parameters().exists("filterN");
		if (filterN) {
			logfile().list("Will filter out sites where reference is 'N'. (argument 'filterN')");
		} else {
			logfile().list("Will keep sites where reference is 'N'. (use 'filterN' to filter out)");
		}
	} else {
		logfile().list("Will identify the most likely among all 6 possible allele combnations. (provide alleles with "
		               "'alleles' or the reference with 'fasta')");
	}

	// write phred-scaled likelihoods?
	const bool usePhredLikelihoods = parameters().exists("phredLik");
	if (usePhredLikelihoods) {
		logfile().list("Will write phred-scaled likelihoods. (parameter phredLik)");
	} else {
		logfile().list("Will write raw likelihoods. (use phredLik to phred-scale)");
	}

	PhredInt minVariantQuality = PhredInt::highest();
	coretools::Probability minMAF{P(0.0)};
	if (parameters().exists("printAll")) {
		logfile().list("Will write all sites and samples. (parameter printAll)");
		if (parameters().exists("minSamplesWithData")) {
			logfile().warning("option 'printAll' overwrites option 'minSamplesWithData' set to ", multiTraverser.minActive(), "!");
		}
		multiTraverser.setMinActive(0);
		minVariantQuality  = PhredInt::highest();
	} else {
		if (multiTraverser.minActive() > 0) {
			logfile().list("Will only print sites for which at least ", multiTraverser.minActive(),
			               " samples have data. (parameter minSamplesWithData)");
		}

		minVariantQuality = parameters().get<PhredInt>("minVariantQual", PhredInt::highest());
		if (minVariantQuality > PhredInt::highest()) {
			logfile().list("Will only print sites with variant quality >= ", minVariantQuality,
			               ". (parameter minVariantQual)");
		}

		minMAF = parameters().get("minMAF", P(0.0));
		if (minMAF > 0.0) {
			logfile().list("Will filter on a minor allele frequency of ", minMAF, ". (parameter 'minMAF')");
		} else {
			logfile().list("Will keep sites regardless of their minor allele frequency. (use 'minMAF' to filter)");
		}
	}

	// limit input
	const size_t limitSites = parameters().get("limitSites", 0);
	logfile().list("Will stop at input position ", limitSites, ". (parameter 'limitSites')");

	// filename tag
	const std::string outname = parameters().get<std::string>("out", "ATLAS_majorMinor");
	logfile().list("Will write output files with tag '" + outname + "'. (parameter 'out')");

	// open vcf file
	const bool writeHeader = !parameters().exists("noVCFHeader");
	if (writeHeader) {
		logfile().list("Will print VCF header. (use 'noVCFHeader' not to print)");
	} else {
		logfile().list("Will not print VCF header. (parameter 'noVCFHeader')");
	}
	genometools::TVCFWriter vcf =
	    coretools::instances::parameters().exists("bgz")
	        ? genometools::TVCFWriter(new GLF::TBGzWriter(outname + ".vcf.gz"), "ATLAS_GLF_Caller",
	                                  multiTraverser.sampleNames(), multiTraverser.chromosomes(), usePhredLikelihoods,
	                                  writeHeader)
	        : genometools::TVCFWriter(outname + ".vcf.gz", "ATLAS_GLF_Caller", multiTraverser.sampleNames(),
	                                  multiTraverser.chromosomes(), usePhredLikelihoods, writeHeader);

	// vars
	logfile().startIndent("Parsing through glf files:");
	coretools::TTimer timer;
	const auto dCounter = std::max<size_t>(1'000'000, 10'000'000 / multiTraverser.numSamples());

	auto &alleles = multiTraverser.alleles();
	genometools::TAlleles::iterator alleleIt;
	for (; !multiTraverser.endOfChrs(); multiTraverser.nextChr()) {
		logfile().startIndent("Parsing chromosome ", multiTraverser.curChr().name(), ":");
		size_t counter  = 0;
		size_t position = dCounter;

		if (alleles) alleleIt = alleles.begin(multiTraverser.curChr().refID());

		for (; !multiTraverser.endOfCurChr(); multiTraverser.nextSite()) {
			const auto &pos = multiTraverser.position();

			if (pos.position() > position) {
				logfile().list("Parsed ", dCounter, " positions in ", timer.formattedTime(), ".");
				position += dCounter;
			}
			MajorMinor::TMMData data;
			if (alleles) { // 1) working with a subset of known alleles
				while (alleleIt != alleles.end() && alleleIt->position < pos) { ++alleleIt; }
				if (alleleIt == alleles.end() || alleleIt->position.refID() > pos.refID()) break;
				DEBUG_ASSERT(alleleIt->position.refID() >= pos.refID());

				if (alleleIt->position.position() > pos.position()) continue;
				DEBUG_ASSERT(alleleIt->position.position() == pos.position());
				data = Estimator::estimate(multiTraverser.site(), maxF, alleleIt->ref, alleleIt->alt, minMAF,
				                           minVariantQuality);
			} else if (fastaReader) { // 2) working with ref
				const auto ref = fastaReader[pos];

				if (filterN && ref == Base::N) continue;
				data = Estimator::estimate(multiTraverser.site(), maxF, ref, minMAF, minVariantQuality);
			} else { // 3) working without external info
				data = Estimator::estimate(multiTraverser.site(), maxF, minMAF, minVariantQuality);
			}
			if (!data.pass) continue;

			vcf.writeSite(multiTraverser.curChr().name(), pos.position(), data.major, data.minor, data.variantQuality,
			              multiTraverser.site());
			++counter;
		}
		const auto tot      = multiTraverser.curChr().to().position();
		const auto filtered = tot - counter;
		logfile().list("Parsed a total of ", tot, " positions, filtered: ", filtered, " (", (100. * filtered) / tot,
					   "%) in ", timer.formattedTime(), ".");
		logfile().endIndent();
	}
	logfile().list("Reached end of glf files in ", timer.formattedTime(), "!");
	logfile().removeIndent();
}

//---------------------------------------------------
// TMajorMinor
//---------------------------------------------------
void TMajorMinor::run() {
	const std::string method = parameters().get<std::string>("method", "MLE");

	const double maxF = parameters().get("maxF", 0.0000001);
	if (method == "Skotte") {
		logfile().list("Will estimate major / minor alleles using the Skotte method with maxF ", maxF,
					   ". (parameters method and maxF)");
		iterate<MajorMinor::TSkotte>(maxF);
	} else if (method == "MLE") {
		logfile().list("Will estimate major / minor alleles using the MLE method with maxF ", maxF,
					   ". (parameters method and maxF)");
		iterate<MajorMinor::TMLE>(maxF);
	} else {
		throw coretools::TUserError("Unknown MajorMinor method '", method, "'!");
	}
}

} // namespace PopulationTools
