//
// Created by madleina on 06.04.22.
//
#include "TSimulator.h"
#include "genometools/VCF/TVcfFile.h"
#include "gtest/gtest.h"

using genometools::Base;

//-------------------------------------------
// TVCFSimulator
//-------------------------------------------

TEST(TVCFSimulator, findMajorMinorAllele) {
	// 1) ref is unique major allele, and there is a unique minor allele
	coretools::TStrongArray<size_t, Base, 4> alleleCounts({10, 20, 5, 15});
	auto ref            = Base::C;
	auto [major, minor] = findMajorMinorAllele(alleleCounts, ref);
	EXPECT_EQ(major, ref);
	EXPECT_EQ(minor, Base::T);

	// 2) ref is unique major allele, and there are multiple minor alleles
	//    -> chose randomly among minor alleles
	alleleCounts = {{0, 0, 10, 0}};
	ref = Base::G;

	size_t N = 1e04;
	coretools::TStrongArray<size_t, Base, 4> counter{};
	for (size_t i = 0; i < N; ++i) {
		auto [major, minor] = findMajorMinorAllele(alleleCounts, ref);
		EXPECT_EQ(major, ref);
		counter[minor]++;
	}
	EXPECT_EQ(counter[Base::A] + counter[Base::C] + counter[Base::T], N);
	EXPECT_TRUE(counter[Base::A] > 0);
	EXPECT_TRUE(counter[Base::C] > 0);
	EXPECT_TRUE(counter[Base::T] > 0);

	// 3) there are multiple major alleles, but ref is none of them
	//    -> this should never happen for bi-allelic simulators!
	alleleCounts = {{10, 0, 10, 0}};
	ref = Base::C;
	EXPECT_ANY_THROW(findMajorMinorAllele(alleleCounts, ref));

	// 4) there are multiple major alleles, one of them is ref
	//    -> chose randomly among major alleles; if ref is not major allele, then it must be the minor allele
	alleleCounts = {{10, 10, 10, 0}};
	ref = Base::G;

	N = 1e04;
	coretools::TStrongArray<size_t, Base, 4> counterMajor{};
	coretools::TStrongArray<size_t, Base, 4> counterMinor{};
	for (size_t i = 0; i < N; ++i) {
		auto [major, minor] = findMajorMinorAllele(alleleCounts, ref);
		EXPECT_TRUE((major == ref || minor == ref));  // must be either major or minor
		EXPECT_FALSE((major == ref && minor == ref)); // can not be both at once
		counterMajor[major]++;
		counterMinor[minor]++;
	}
	EXPECT_EQ(counterMajor[Base::A] + counterMajor[Base::C] +
	              counterMajor[Base::G],
	          N);
	EXPECT_EQ(counterMinor[Base::A] + counterMinor[Base::C] +
	              counterMinor[Base::G],
	          N);

	EXPECT_TRUE(counterMajor[Base::A] > 0);
	EXPECT_TRUE(counterMajor[Base::G] > 0);
	EXPECT_TRUE(counterMajor[Base::C] > 0);
	EXPECT_TRUE(counterMinor[Base::A] > 0);
	EXPECT_TRUE(counterMinor[Base::G] > 0);
	EXPECT_TRUE(counterMinor[Base::C] > 0);
}

void simulate(){
	Simulations::TVCFSimulator vcfSimulator("HW");
	vcfSimulator.runSimulations();
}

TEST(TVCFSimulator, integrationTest) {
	coretools::instances::randomGenerator().setSeed(0);
	size_t numSamples = 20;
	coretools::instances::parameters().add("sampleSize", numSamples);
	coretools::instances::parameters().add("chrLength", "1000,1000");
	coretools::instances::parameters().add("ploidy", "1,2");
	coretools::instances::parameters().add("depth", "10");
	coretools::instances::parameters().add("refN", "0");
	coretools::instances::parameters().add("out", "./vcfSimulatorTest");

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
