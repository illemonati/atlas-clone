//
// Created by caduffm on 10/27/21.
//
#include "../TVcfConverter.h"
#include "gtest/gtest.h"

#include <cstdint>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "coretools/Files/TFile.h"
#include "coretools/Files/gzstream.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Types/probability.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "genometools/TSampleLikelihoods.h"

using namespace testing;
using namespace coretools::instances;
using namespace coretools::str;
using namespace coretools;
using namespace genometools;

class TVCFConverterTest : public Test {
protected:
	std::vector<uint16_t> phred_g1         = {0, 1, 0, 0, 3, 0, 1, 2, 0, 0, 2, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 3, 1,
                                      0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 2, 3, 2, 0, 3, 0, 0, 1, 0, 0, 0, 1, 2, 0};
	std::vector<uint16_t> phred_g2         = {1, 2, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 2, 1, 2, 2, 0, 1, 1, 0, 2, 0, 0, 0,
                                      0, 0, 0, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1};
	std::vector<uint16_t> phred_g3         = {1, 0, 1, 3, 1, 2, 0, 1, 0, 1, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                                      1, 3, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 2, 3, 0, 0, 1, 0, 2, 1, 0, 1};
	std::vector<uint16_t> depth            = {14, 11, 9, 9, 6,  8, 12, 9,  10, 11, 5,  14, 7,  13, 13, 11, 11,
                                   5,  9,  8, 9, 8,  6, 7,  10, 14, 10, 11, 13, 10, 7,  9,  11, 13,
                                   8,  7,  9, 6, 10, 6, 15, 15, 10, 9,  11, 9,  6,  14, 9,  14};
	std::vector<std::string> sampleNames   = {"Indiv1", "Indiv2", "Indiv3", "Indiv4", "Indiv5"};
	std::vector<std::string> samplesToKeep = {"Indiv4", "Indiv1", "Indiv5"}; // omit Indiv3 and Indiv4; and shuffle
	std::vector<size_t> indexInSampleNames = {3, 0, 4};

	static void _writeHeader(gz::ogzstream &VCF, const std::vector<std::string> &SampleNames) {
		// write info
		VCF << "##fileformat=VCFv4.2\n";
		VCF << "##source=Simulation\n";
		VCF << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		VCF << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
		VCF << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">\n";
		VCF << "##FORMAT=<ID=DP,Number=G,Type=Integer,Description=\"Depth at site\">\n";

		// write header
		VCF << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
		for (const auto &SampleName : SampleNames) { VCF << "\t" << SampleName; }
	}

	static void _writeCell(gz::ogzstream &VCF, const TSampleLikelihoods<PhredIntProbability> &SampleLikelihoods,
	                       size_t Depth) {
		using BG = BiallelicGenotype;

		// genotype with lowest phred-score is the observed genotype
		BG observedGenotype   = SampleLikelihoods.mostLikelyGenotype();
		BG secondBestGenotype = SampleLikelihoods.secondMostLikelyGenotype();

		// write to vcf
		VCF << "\t" << toString(observedGenotype) << ":" << SampleLikelihoods[secondBestGenotype] << ":";
		VCF << toString(SampleLikelihoods[BG::homoFirst]) + "," + toString(SampleLikelihoods[BG::het]) + "," +
		           toString(SampleLikelihoods[BG::homoSecond]);
		VCF << ":" << Depth;
	};

public:
	std::string filename;

	const uint32_t numLoci  = 10;
	const uint32_t numIndiv = 5;

	void SetUp() override {
		parameters().clear();

		filename = "test.vcf.gz";
		parameters().addParameter("vcf", filename);
		parameters().addParameter("minMAF", "0");
	}

	void writeVcfFile() {
		gz::ogzstream vcf;
		vcf.open(filename.c_str());
		_writeHeader(vcf, sampleNames);

		std::string chr    = "junk1";
		size_t linearIndex = 0;
		for (size_t l = 0; l < numLoci; ++l) {
			vcf << '\n' << chr << '\t' << l + 1 << "\t.\tA\tC\t50\t.\t.\tGT:GQ:PL:DP";
			for (size_t i = 0; i < numIndiv; ++i, ++linearIndex) {
				// write to vcf
				auto GTL0 = PhredIntProbability(phred_g1[linearIndex]);
				auto GTL1 = PhredIntProbability(phred_g2[linearIndex]);
				auto GTL2 = PhredIntProbability(phred_g3[linearIndex]);
				_writeCell(vcf, {GTL0, GTL1, GTL2}, depth[linearIndex]);
			}
		}
		vcf << "\n";
		vcf.close();
	}

	void writeSampleFile() {
		TOutputFile sampleFile("test.samples");
		for (auto &it : samplesToKeep) { sampleFile.writeln(it); }
		// add to parameters
		parameters().addParameter("samples", "test.samples");
	}
};

TEST_F(TVCFConverterTest, beagle) {
	writeVcfFile();

	// convert to beagle
	VCF::TVcfToBeagle converter;
	converter.run();

	// read beagle
	TInputFile beagle("test.beagle.gz", TFile_Filetype::header);
	EXPECT_EQ(beagle.headerAt(0), "marker");
	EXPECT_EQ(beagle.headerAt(1), "allele1");
	EXPECT_EQ(beagle.headerAt(2), "allele2");
	for (size_t i = 0; i < numIndiv; ++i) {
		EXPECT_EQ(beagle.headerAt(3 + i * 3), sampleNames[i]);
		EXPECT_EQ(beagle.headerAt(3 + i * 3 + 1), sampleNames[i]);
		EXPECT_EQ(beagle.headerAt(3 + i * 3 + 2), sampleNames[i]);
	}

	// check if genotype likelihoods are as expected
	std::vector<std::string> line;
	size_t linearIndex = 0;
	for (size_t l = 0; l < numLoci; ++l) {
		beagle.read(line);
		EXPECT_EQ(line[0], "junk1_" + toString(l + 1));
		EXPECT_EQ(line[1], "A");
		EXPECT_EQ(line[2], "C");
		for (size_t i = 0; i < numIndiv; ++i, ++linearIndex) {
			const double sum = static_cast<Probability>(PhredIntProbability(phred_g1[linearIndex])).get() +
			                   (Probability)PhredIntProbability(phred_g2[linearIndex]) +
			                   (Probability)PhredIntProbability(phred_g3[linearIndex]);
			const double gtl1 = fromString<Probability>(line[3 + i * 3]);
			const double gtl2 = fromString<Probability>(line[3 + i * 3 + 1]);
			const double gtl3 = fromString<Probability>(line[3 + i * 3 + 2]);

			// genotype likelihoods (we loose some precision when reading/writing, so only expect near)
			EXPECT_NEAR(gtl1, (Probability)PhredIntProbability(phred_g1[linearIndex]) / sum, 0.00001);
			EXPECT_NEAR(gtl2, (Probability)PhredIntProbability(phred_g2[linearIndex]) / sum, 0.00001);
			EXPECT_NEAR(gtl3, (Probability)PhredIntProbability(phred_g3[linearIndex]) / sum, 0.00001);
			EXPECT_NEAR(gtl1 + gtl2 + gtl3, 1.0, 0.00001);
		}
	}
}

TEST_F(TVCFConverterTest, beagle_withSamples) {
	writeVcfFile();
	writeSampleFile();

	// convert to beagle
	VCF::TVcfToBeagle converter;
	converter.run();

	// read beagle
	TInputFile beagle("test.beagle.gz", TFile_Filetype::header);
	EXPECT_EQ(beagle.headerAt(0), "marker");
	EXPECT_EQ(beagle.headerAt(1), "allele1");
	EXPECT_EQ(beagle.headerAt(2), "allele2");
	for (size_t i = 0; i < samplesToKeep.size(); ++i) {
		EXPECT_EQ(beagle.headerAt(3 + i * 3), samplesToKeep[i]);
		EXPECT_EQ(beagle.headerAt(3 + i * 3 + 1), samplesToKeep[i]);
		EXPECT_EQ(beagle.headerAt(3 + i * 3 + 2), samplesToKeep[i]);
	}

	// check if genotype likelihoods are as expected
	std::vector<std::string> line;
	for (size_t l = 0; l < numLoci; ++l) {
		beagle.read(line);
		EXPECT_EQ(line[0], "junk1_" + toString(l + 1));
		EXPECT_EQ(line[1], "A");
		EXPECT_EQ(line[2], "C");
		for (size_t i = 0; i < samplesToKeep.size(); ++i) {
			size_t relevantIndex = l * numIndiv + indexInSampleNames[i];
			// genotype likelihoods (we loose some precision when reading/writing, so only expect near)
			const double sum     = ((Probability)PhredIntProbability(phred_g1[relevantIndex])).get() +
			                   (Probability)PhredIntProbability(phred_g2[relevantIndex]) +
			                   (Probability)PhredIntProbability(phred_g3[relevantIndex]);
			const double gtl1 = fromString<Probability>(line[3 + i * 3]);
			const double gtl2 = fromString<Probability>(line[3 + i * 3 + 1]);
			const double gtl3 = fromString<Probability>(line[3 + i * 3 + 2]);
			EXPECT_NEAR(gtl1, (Probability)PhredIntProbability(phred_g1[relevantIndex]) / sum, 0.00001);
			EXPECT_NEAR(gtl2, (Probability)PhredIntProbability(phred_g2[relevantIndex]) / sum, 0.00001);
			EXPECT_NEAR(gtl3, (Probability)PhredIntProbability(phred_g3[relevantIndex]) / sum, 0.00001);
			EXPECT_NEAR(gtl1 + gtl2 + gtl3, 1.0, 0.00001);
		}
	}
}

TEST_F(TVCFConverterTest, geno) {
	writeVcfFile();

	// convert to geno
	VCF::TVcfToGeno converter;
	converter.run();

	// read geno
	TInputFile geno("test.geno", TFile_Filetype::fixed);

	// check if genotypes are as expected
	std::vector<std::string> line;
	for (size_t l = 0; l < numLoci; ++l) {
		geno.read(line);
		std::vector<char> genotypes(line[0].begin(), line[0].end());
		for (size_t i = 0; i < numIndiv; ++i) {
			// genotypes
			size_t relevantIndex = l * numIndiv + i;
			auto GTL0            = PhredIntProbability(phred_g1[relevantIndex]);
			auto GTL1            = PhredIntProbability(phred_g2[relevantIndex]);
			auto GTL2            = PhredIntProbability(phred_g3[relevantIndex]);
			TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
			BiallelicGenotype observedGenotype = sampleLikelihoods.mostLikelyGenotype();

			EXPECT_EQ(fromString<uint8_t>(std::string(1, genotypes[i])),
			          2 - (uint8_t)observedGenotype); // genotype defined based on reference allele
		}
	}
}

TEST_F(TVCFConverterTest, geno_withSamples) {
	writeVcfFile();
	writeSampleFile();

	// convert to geno
	VCF::TVcfToGeno converter;
	converter.run();

	// read geno
	TInputFile geno("test.geno", TFile_Filetype::fixed);

	// check if genotypes are as expected
	std::vector<std::string> line;
	for (size_t l = 0; l < numLoci; ++l) {
		geno.read(line);
		std::vector<char> genotypes(line[0].begin(), line[0].end());
		for (size_t i = 0; i < samplesToKeep.size(); ++i) {
			// genotypes
			size_t relevantIndex = l * numIndiv + indexInSampleNames[i];
			auto GTL0            = PhredIntProbability(phred_g1[relevantIndex]);
			auto GTL1            = PhredIntProbability(phred_g2[relevantIndex]);
			auto GTL2            = PhredIntProbability(phred_g3[relevantIndex]);
			TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
			BiallelicGenotype observedGenotype = sampleLikelihoods.mostLikelyGenotype();

			// genotype defined based on reference allele
			EXPECT_EQ(fromString<uint8_t>(std::string(1, genotypes[i])), 2 - (uint8_t)observedGenotype);
		}
	}
}

TEST_F(TVCFConverterTest, lfmmCalledGeno) {
	writeVcfFile();

	// convert to lfmm
	VCF::TVcfToLFMM<true> converter;
	converter.run();

	// read lfmm
	TInputFile lfmm("test.lfmm", TFile_Filetype::fixed);

	// check if genotypes are as expected
	std::vector<std::string> line;
	for (size_t i = 0; i < numIndiv; ++i) {
		lfmm.read(line);
		for (size_t l = 0; l < numLoci; ++l) {
			// genotypes
			size_t relevantIndex = l * numIndiv + i;
			auto GTL0            = PhredIntProbability(phred_g1[relevantIndex]);
			auto GTL1            = PhredIntProbability(phred_g2[relevantIndex]);
			auto GTL2            = PhredIntProbability(phred_g3[relevantIndex]);
			TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
			BiallelicGenotype observedGenotype = sampleLikelihoods.mostLikelyGenotype();

			EXPECT_EQ(fromString<uint8_t>(line[l]), (uint8_t)observedGenotype);
		}
	}
}

TEST_F(TVCFConverterTest, lfmmCalledGeno_withSamples) {
	writeVcfFile();
	writeSampleFile();

	// convert to lfmm
	VCF::TVcfToLFMM<true> converter;
	converter.run();

	// read lfmm
	TInputFile lfmm("test.lfmm", TFile_Filetype::fixed);

	// check if genotypes are as expected
	std::vector<std::string> line;
	for (size_t i = 0; i < samplesToKeep.size(); ++i) {
		lfmm.read(line);
		for (size_t l = 0; l < numLoci; ++l) {
			// genotypes
			size_t relevantIndex = l * numIndiv + indexInSampleNames[i];
			auto GTL0            = PhredIntProbability(phred_g1[relevantIndex]);
			auto GTL1            = PhredIntProbability(phred_g2[relevantIndex]);
			auto GTL2            = PhredIntProbability(phred_g3[relevantIndex]);
			TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
			BiallelicGenotype observedGenotype = sampleLikelihoods.mostLikelyGenotype();

			EXPECT_EQ(fromString<uint8_t>(line[l]), (uint8_t)observedGenotype);
		}
	}
}

TEST_F(TVCFConverterTest, lfmmMeanPosteriorGeno) {
	writeVcfFile();

	// convert to lfmm
	VCF::TVcfToLFMM<false> converter;
	converter.run();

	// read lfmm
	TInputFile lfmm("test.lfmm", TFile_Filetype::fixed);

	// check if genotype likelihoods are as expected
	std::vector<std::string> line;
	for (size_t i = 0; i < numIndiv; ++i) {
		lfmm.read(line);
		for (size_t l = 0; l < numLoci; ++l) {
			// genotypes
			size_t relevantIndex = l * numIndiv + i;
			auto GTL0            = PhredIntProbability(phred_g1[relevantIndex]);
			auto GTL1            = PhredIntProbability(phred_g2[relevantIndex]);
			auto GTL2            = PhredIntProbability(phred_g3[relevantIndex]);
			TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
			double posteriorGenotype = sampleLikelihoods.meanPosteriorGenotype();

			EXPECT_NEAR(fromString<double>(line[l]), posteriorGenotype, 0.00001);
		}
	}
}

TEST_F(TVCFConverterTest, lfmmMeanPosteriorGeno_withSamples) {
	writeVcfFile();
	writeSampleFile();

	// convert to lfmm
	VCF::TVcfToLFMM<false> converter;
	converter.run();

	// read lfmm
	TInputFile lfmm("test.lfmm", TFile_Filetype::fixed);

	// check if genotypes are as expected
	std::vector<std::string> line;
	for (size_t i = 0; i < samplesToKeep.size(); ++i) {
		lfmm.read(line);
		for (size_t l = 0; l < numLoci; ++l) {
			// genotypes
			size_t relevantIndex = l * numIndiv + indexInSampleNames[i];
			auto GTL0            = PhredIntProbability(phred_g1[relevantIndex]);
			auto GTL1            = PhredIntProbability(phred_g2[relevantIndex]);
			auto GTL2            = PhredIntProbability(phred_g3[relevantIndex]);
			TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
			double posteriorGenotype = sampleLikelihoods.meanPosteriorGenotype();

			EXPECT_NEAR(fromString<double>(line[l]), posteriorGenotype, 0.00001);
		}
	}
}

TEST_F(TVCFConverterTest, sambada) {
	writeVcfFile();

	// convert to lfmm
	VCF::TVcfToSambada converter;
	converter.run();

	// read lfmm
	TInputFile sambada("test.sambada", TFile_Filetype::header);

	// check if genotypes are as expected
	std::vector<std::string> line;
	for (size_t i = 0; i < numIndiv; ++i) {
		sambada.read(line);
		EXPECT_EQ(line[0], sampleNames[i]);
		for (size_t l = 0; l < numLoci; ++l) {
			// genotypes
			size_t relevantIndex = l * numIndiv + i;
			auto GTL0            = PhredIntProbability(phred_g1[relevantIndex]);
			auto GTL1            = PhredIntProbability(phred_g2[relevantIndex]);
			auto GTL2            = PhredIntProbability(phred_g3[relevantIndex]);
			TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
			BiallelicGenotype observedGenotype = sampleLikelihoods.mostLikelyGenotype();

			for (size_t g = 0; g < 3; ++g) {
				if (g == (uint8_t)observedGenotype) {
					EXPECT_EQ(fromString<size_t>(line[3 * l + g + 1]), 1);
				} else {
					EXPECT_EQ(fromString<size_t>(line[3 * l + g + 1]), 0);
				}
			}
		}
	}
}

TEST_F(TVCFConverterTest, sambada_withSamples) {
	writeVcfFile();
	writeSampleFile();

	// convert to lfmm
	VCF::TVcfToSambada converter;
	converter.run();

	// read sambada
	TInputFile sambada("test.sambada", TFile_Filetype::header);

	// check if genotypes are as expected
	std::vector<std::string> line;
	for (size_t i = 0; i < samplesToKeep.size(); ++i) {
		sambada.read(line);
		EXPECT_EQ(line[0], samplesToKeep[i]);
		for (size_t l = 0; l < numLoci; ++l) {
			// genotypes
			size_t relevantIndex = l * numIndiv + indexInSampleNames[i];
			auto GTL0            = PhredIntProbability(phred_g1[relevantIndex]);
			auto GTL1            = PhredIntProbability(phred_g2[relevantIndex]);
			auto GTL2            = PhredIntProbability(phred_g3[relevantIndex]);
			TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
			BiallelicGenotype observedGenotype = sampleLikelihoods.mostLikelyGenotype();

			for (size_t g = 0; g < 3; ++g) {
				if (g == (uint8_t)observedGenotype) {
					EXPECT_EQ(fromString<uint8_t>(line[3 * l + g + 1]), 1);
				} else {
					EXPECT_EQ(fromString<uint8_t>(line[3 * l + g + 1]), 0);
				}
			}
		}
	}
}

TEST_F(TVCFConverterTest, vcfToPosFile) {
	writeVcfFile();
	writeSampleFile();

	// convert to lfmm
	VCF::TVcfToPosFile converter;
	converter.run();

	// read pos file
	TInputFile pos("test.pos", TFile_Filetype::fixed);

	// check if positions are as expected
	std::vector<std::string> line;
	for (size_t l = 0; l < numLoci; ++l) {
		pos.read(line);
		EXPECT_EQ(line[0], "junk1");
		EXPECT_EQ(line[1], toString(l + 1));
		EXPECT_EQ(line[2], "A");
		EXPECT_EQ(line[3], "C");
	}
}
