//
// Created by madleina on 06.04.22.
//
#include "../TSimulator.h"
#include "gtest/gtest.h"

#include <stddef.h>
#include <istream>
#include <stdexcept>
#include <string>

#include "../genometools/core/VCF/TVcfFile.h"
#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "SFS.h"
#include "TGlfMultiReader.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TStrongArray.h"

using namespace testing;
using namespace genometools;
using namespace Simulations;

//-------------------------------------------
// TVCFSimulator
//-------------------------------------------

class TVCFSimulatorBridge : public TVCFSimulator {
	// make protected functions public
public:
	TVCFSimulatorBridge() : TVCFSimulator("HW") {}

	GLF::TMultiGLFDataSampleOneAllelicCombination _calculateGenotypeLikelihoods(size_t NumRef, size_t NumAlt,
	                                                                            bool IsDiploid) {
		return TVCFSimulator::_calculateGenotypeLikelihoods(NumRef, NumAlt, IsDiploid);
	}
	auto _findMajorMinorAllele(coretools::TStrongArray<size_t, genometools::Base, 4> AlleleCounts,
	                           genometools::Base RefAllele) {
		return TVCFSimulator::_findMajorMinorAllele(AlleleCounts, RefAllele);
	}
};

TEST(TVCFSimulator, findMajorMinorAllele_1) {
	TVCFSimulatorBridge vcfSimulator;

	// 1) ref is unique major allele, and there is a unique minor allele
	coretools::TStrongArray<size_t, genometools::Base, 4> alleleCounts({10, 20, 5, 15});
	auto ref            = genometools::Base::C;
	auto [major, minor] = vcfSimulator._findMajorMinorAllele(alleleCounts, ref);
	EXPECT_EQ(major, ref);
	EXPECT_EQ(minor, genometools::Base::T);
}

TEST(TVCFSimulator, findMajorMinorAllele_2) {
	TVCFSimulatorBridge vcfSimulator;

	// 2) ref is unique major allele, and there are multiple minor alleles
	//    -> chose randomly among minor alleles
	coretools::TStrongArray<size_t, genometools::Base, 4> alleleCounts({0, 0, 10, 0});
	auto ref = genometools::Base::G;

	size_t N = 1e04;
	coretools::TStrongArray<size_t, genometools::Base, 4> counter{};
	for (size_t i = 0; i < N; ++i) {
		auto [major, minor] = vcfSimulator._findMajorMinorAllele(alleleCounts, ref);
		EXPECT_EQ(major, ref);
		counter[minor]++;
	}
	EXPECT_EQ(counter[genometools::Base::A] + counter[genometools::Base::C] + counter[genometools::Base::T], N);
	EXPECT_TRUE(counter[genometools::Base::A] > 0);
	EXPECT_TRUE(counter[genometools::Base::C] > 0);
	EXPECT_TRUE(counter[genometools::Base::T] > 0);
}

TEST(TVCFSimulator, findMajorMinorAllele_3) {
	TVCFSimulatorBridge vcfSimulator;

	// 3) there are multiple major alleles, but ref is none of them
	//    -> this should never happen for bi-allelic simulators!
	coretools::TStrongArray<size_t, genometools::Base, 4> alleleCounts({10, 0, 10, 0});
	auto ref = genometools::Base::C;
	EXPECT_THROW(vcfSimulator._findMajorMinorAllele(alleleCounts, ref), std::runtime_error);
}

TEST(TVCFSimulator, findMajorMinorAllele_4) {
	TVCFSimulatorBridge vcfSimulator;

	// 4) there are multiple major alleles, one of them is ref
	//    -> chose randomly among major alleles; if ref is not major allele, then it must be the minor allele
	coretools::TStrongArray<size_t, genometools::Base, 4> alleleCounts({10, 10, 10, 0});
	auto ref = genometools::Base::G;

	size_t N = 1e04;
	coretools::TStrongArray<size_t, genometools::Base, 4> counterMajor{};
	coretools::TStrongArray<size_t, genometools::Base, 4> counterMinor{};
	for (size_t i = 0; i < N; ++i) {
		auto [major, minor] = vcfSimulator._findMajorMinorAllele(alleleCounts, ref);
		EXPECT_TRUE((major == ref || minor == ref));  // must be either major or minor
		EXPECT_FALSE((major == ref && minor == ref)); // can not be both at once
		counterMajor[major]++;
		counterMinor[minor]++;
	}
	EXPECT_EQ(counterMajor[genometools::Base::A] + counterMajor[genometools::Base::C] +
	              counterMajor[genometools::Base::G],
	          N);
	EXPECT_EQ(counterMinor[genometools::Base::A] + counterMinor[genometools::Base::C] +
	              counterMinor[genometools::Base::G],
	          N);

	EXPECT_TRUE(counterMajor[genometools::Base::A] > 0);
	EXPECT_TRUE(counterMajor[genometools::Base::G] > 0);
	EXPECT_TRUE(counterMajor[genometools::Base::C] > 0);
	EXPECT_TRUE(counterMinor[genometools::Base::A] > 0);
	EXPECT_TRUE(counterMinor[genometools::Base::G] > 0);
	EXPECT_TRUE(counterMinor[genometools::Base::C] > 0);
}

TEST(TVCFSimulator, _calculateGenotypeLikelihoods_diploid) {
	coretools::instances::parameters().addParameter("error", 0.05);
	TVCFSimulatorBridge vcfSimulator;

	auto GTL = vcfSimulator._calculateGenotypeLikelihoods(10, 20, true);
	EXPECT_EQ(HighPrecisionPhredIntProbability(GTL[genometools::BiallelicGenotype::homoFirst]), 26243);
	EXPECT_EQ(HighPrecisionPhredIntProbability(GTL[genometools::BiallelicGenotype::het]), 9031);
	EXPECT_EQ(HighPrecisionPhredIntProbability(GTL[genometools::BiallelicGenotype::homoSecond]), 13456);
}

TEST(TVCFSimulator, _calculateGenotypeLikelihoods_haploid) {
	coretools::instances::parameters().addParameter("error", 0.05);
	TVCFSimulatorBridge vcfSimulator;

	auto GTL = vcfSimulator._calculateGenotypeLikelihoods(10, 20, false);
	EXPECT_EQ(HighPrecisionPhredIntProbability(GTL[genometools::BiallelicGenotype::haploidFirst]), 26243);
	EXPECT_EQ(HighPrecisionPhredIntProbability(GTL[genometools::BiallelicGenotype::haploidSecond]), 13456);
}

void simulate(){
	TVCFSimulator vcfSimulator("HW");
	vcfSimulator.runSimulations();
}

TEST(TVCFSimulator, integrationTest) {
	coretools::instances::randomGenerator().setSeed(0);
	size_t numSamples = 20;
	coretools::instances::parameters().addParameter("sampleSize", numSamples);
	coretools::instances::parameters().addParameter("chrLength", "1000,1000");
	coretools::instances::parameters().addParameter("ploidy", "1,2");
	coretools::instances::parameters().addParameter("out", "./vcfSimulatorTest");

	// simulate
	simulate();

	// read reference
	std::string reference, tmp;
	std::ifstream fastaFile("./vcfSimulatorTest.fasta");
	while (std::getline(fastaFile, tmp)) {
		if (tmp.empty()) { continue; }
		if (tmp[0] != '>') { reference += tmp; }
	}

	// read
	genometools::TVcfFileSingleLine vcfFile;
	vcfFile.openStream("./vcfSimulatorTest.vcf.gz");

	// enable parsers
	vcfFile.enablePositionParsing();
	vcfFile.enableVariantParsing();
	vcfFile.enableVariantQualityParsing();
	vcfFile.enableFormatParsing();
	vcfFile.enableSampleParsing();

	EXPECT_EQ(vcfFile.numSamples(), numSamples);

	size_t c     = 0;
	size_t refID = 0;
	size_t pos   = 0;
	while (vcfFile.next()) {
		if (c == 1000) { // second chromosome
			pos   = 0;
			refID = 1;
		}

		// check for reference allele
		EXPECT_EQ(vcfFile.getRefAllele(), std::string(1, reference[c]));
		EXPECT_TRUE(vcfFile.getRefAllele() != vcfFile.getFirstAltAllele());
		EXPECT_TRUE(vcfFile.getNumAlleles() == 2);
		EXPECT_EQ(vcfFile.getNumAlleles(), 2); // always biallelic
		EXPECT_TRUE(vcfFile.isBialleleicSNP());

		// check chromosome
		if (refID == 0) {
			EXPECT_EQ(vcfFile.chr(), "chr1");
			EXPECT_EQ(vcfFile.position(), pos + 1);
		} else {
			EXPECT_EQ(vcfFile.chr(), "chr2");
			EXPECT_EQ(vcfFile.position(), pos + 1);
		}

		for (size_t i = 0; i < numSamples; i++) {
			if (vcfFile.sampleIsMissing(i)) { continue; }

			// check depth: never zero (then it should be missing)
			EXPECT_TRUE(vcfFile.sampleDepth(i) > 0);

			// check haploid/diploid
			if (refID == 0) {
				EXPECT_TRUE(vcfFile.sampleIsHaploid(i));
				EXPECT_FALSE(vcfFile.sampleIsDiploid(i));
			} else {
				EXPECT_FALSE(vcfFile.sampleIsHaploid(i));
				EXPECT_TRUE(vcfFile.sampleIsDiploid(i));
			}
		}

		++c;
		++pos;
	}
}
