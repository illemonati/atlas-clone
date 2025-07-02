/*
 * TMajorMinor.cpp
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#include "TMajorMinor.h"
#include "TSkotte.h"
#include "TMajorMinorMLE.h"

#include <algorithm>


#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Types/probability.h"
#include "genometools/GLF/GLF.h"
#include "genometools/Genotypes/TFrequencies.h"
#include "genometools/VCF/TVCFWriter.h"

#include "genometools/TAlleles.h"
#include "TBgzWriter.h"
#include "genometools/GLF/TGLFMultiReader.h"

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
} // namespace MajorMinor

template<typename Estimator> void iterate(double maxF) {
	// open GLF files
	genometools::TGLFMultiReader glfReader;

	// use known alleles or reference allele, if provided

	bool filterN = false;

	genometools::TAlleles alleles;
	if (parameters().exists("alleles")) {
		logfile().startIndent("Will limit analysis to sites with known alleles (parameter 'alleles'):");
		const auto filename = parameters().get("alleles");
		alleles.parse(filename, glfReader.chromosomes());
		glfReader.setAlleles(alleles);
		logfile().endIndent();
	} else if (parameters().exists("fasta")) {
		logfile().list("Will use reference allele and only identify the most likely alternative allele. (argument: fasta)");
		const std::string fastaFile = parameters().get<std::string>("fasta");
		logfile().list("Reading reference sequence from '" + fastaFile + "'");
		glfReader.addReference(fastaFile);
		filterN = parameters().exists("filterN");
		if (filterN) {
			logfile().list("Will filter out sites where reference is 'N'. (argument 'filterN')");
		} else {
			logfile().list("Will keep sites where reference is 'N'. (use 'filterN' to filter out)");
		}
	} else {
		logfile().list("Will identify the most likely among all 6 possible allele combnations. (provide alleles with 'alleles' or the reference with 'fasta')");
	}

	// write phred-scaled likelihoods?
	const bool usePhredLikelihoods = parameters().exists("phredLik");
	if (usePhredLikelihoods) {
		logfile().list("Will write phred-scaled likelihoods. (parameter phredLik)");
	} else {
		logfile().list("Will write raw likelihoods. (use phredLik to phred-scale)");
	}

	// read filters

	const auto hasRef  = glfReader.hasRef();
	if (filterN) {
		logfile().list("Will filter out sites where reference is 'N'. (argument filterN)");
	}

	size_t minSamplesWithData = 1;
	PhredInt minVariantQuality = PhredInt::highest();
	coretools::Probability minMAF{P(0.0)};
	if (parameters().exists("printAll")) {
		logfile().list("Will write all sites and samples. (parameter printAll)");
		minSamplesWithData = 0;
		minVariantQuality  = PhredInt::highest();
	} else {
		minSamplesWithData = parameters().get<size_t>("minSamplesWithData", 1);
		if (minSamplesWithData > 0) {
			logfile().list("Will only print sites for which at least ", minSamplesWithData,
						   " samples have data. (parameter minSamplesWithData)");
		}

		minVariantQuality = parameters().get<PhredInt>(
			"minVariantQual", PhredInt::highest());
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
	glfReader.setMinSamplesWithData(minSamplesWithData);

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
									  glfReader.sampleNames(), glfReader.chromosomes(), usePhredLikelihoods, writeHeader)
			: genometools::TVCFWriter(outname + ".vcf.gz", "ATLAS_GLF_Caller", glfReader.sampleNames(),
									  glfReader.chromosomes(), usePhredLikelihoods, writeHeader);

	// vars
	logfile().startIndent("Parsing through glf files:");
	coretools::TTimer timer;
	constexpr size_t dCounter = 1000000;
	size_t counter            = 0;
	size_t counterF           = 0;
	size_t nextPrint          = dCounter;

	std::vector<MajorMinor::TMMData> data;
	for (; !glfReader.empty(); glfReader.popFront()) {
		const auto& ids = glfReader.front();
		data.assign(ids.back() + 1, MajorMinor::failedTMMData);

		if (alleles) {
			// 1) when working with a subset of known alleles
			const auto begin = alleles.begin(glfReader.curWindow());
			size_t iId = 0;
			for (auto it = begin; it != alleles.end() && it->position < glfReader.curWindow().to(); ++it) {
				const auto iW = it->position.position() - glfReader.curWindow().from().position();
				while (iId < ids.size() && ids[iId] < iW) ++iId;
				if (iId == ids.size()) break;
				if (ids[iId] == iW) {
					data[iW] = Estimator::estimate(glfReader.data(iW), maxF, it->ref, it->alt, minMAF, minVariantQuality);
				}
			}
		} else if (hasRef) {
			// 2) when working with ref
			const auto refs = glfReader.refView();
			for (size_t i = 0; i < ids.size(); ++i) {
				const auto iW = ids[i];
				if (filterN && refs[iW] == Base::N) {
					data[iW].pass = false;
				} else {
					data[iW] = Estimator::estimate(glfReader.data(iW), maxF, refs[iW], minMAF, minVariantQuality);
				}
			}
		} else {
			// 3) working with raw data / no external info
			for (size_t i = 0; i < ids.size(); ++i) {
				const auto iW = ids[i];
				data[iW] = Estimator::estimate(glfReader.data(iW), maxF, minMAF, minVariantQuality);
			}
		}

		// pass filter?
		for (size_t i = 0; i < ids.size(); ++i) {
			const auto iW  = ids[i];
			const auto &di = data[iW];
			if (!di.pass) {
				++counterF;
				continue;
			}

			// write to VCF
			vcf.writeSite(glfReader.curChrName(), glfReader.position(iW).position(), di.major,
						  di.minor, di.variantQuality, glfReader.data(iW));
		}
		counter  += ids.size();
		counterF += (glfReader.curWindow().size() - ids.size());

		// report progress
		if (counter >= nextPrint) {
			logfile().list("Parsed ", nextPrint, " positions in ", timer.formattedTime(), ".");
			while (nextPrint <= counter) nextPrint += dCounter;
		}

		if (limitSites > 0 && counter == limitSites) break;
	}

	logfile().list("Reached end of glf files!");
	logfile().list("Parsed a total of ", counter, " positions, filtered: ", counterF, " (", (100.*counterF)/counter, "%).");
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
