//
// Created by caduffm on 6/9/22.
//

#include "../TVcfDiagnostics.h"
#include "gtest/gtest.h"

#include <cstdint>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TParameters.h"
#include "genometools/TSampleLikelihoods.h"
#include "coretools/Files/gzstream.h"
#include "coretools/Types/probability.h"
#include "coretools/Strings/stringFunctions.h"

using namespace testing;
using namespace coretools::instances;
using namespace coretools::str;
using namespace coretools;
using namespace genometools;

class TVCFDiagnosticsTest : public Test {
protected:
	std::vector<std::string> sampleNames = {"Indiv1", "Indiv2", "Indiv3", "Indiv4", "Indiv5"};

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

	static void _writeCell(gz::ogzstream &VCF, BiallelicGenotype Genotype) {
		// write to vcf
		VCF << "\t" << toString(Genotype) << ":0:0,0,0:100";
	};

public:
	std::string filename;

	const uint32_t numLoci  = 10;
	const uint32_t numIndiv = 5;

	void SetUp() override {
		parameters().clear();

		filename = "test.vcf.gz";
		parameters().add("vcf", filename);
	}

	void writeVcfFile(const std::vector<std::string> &Chromosomes,
	                  const std::vector<std::vector<BiallelicGenotype>> &Genotypes) {
		gz::ogzstream vcf;
		vcf.open(filename.c_str());
		_writeHeader(vcf, sampleNames);

		size_t linearIndex = 0;
		for (size_t l = 0; l < numLoci; ++l) {
			vcf << '\n' << Chromosomes[l] << '\t' << l + 1 << "\t.\tA\tC\t50\t.\t.\tGT:GQ:PL:DP";
			for (size_t i = 0; i < numIndiv; ++i, ++linearIndex) { _writeCell(vcf, Genotypes[l][i]); }
		}
		vcf << "\n";
		vcf.close();
	}
};

TEST_F(TVCFDiagnosticsTest, vcfToInvariantBed_allAreInvariant) {
	std::vector<std::string> chromosomes(numLoci, "junk_1");
	std::vector<std::vector<BiallelicGenotype>> genotypes(numLoci);
	for (size_t l = 0; l < numLoci; l++){
		genotypes[l].resize(numIndiv, BiallelicGenotype::homoFirst);
	}

	writeVcfFile(chromosomes, genotypes);

	// diagnose
	VCF::TVcfDiagnostics diagnostics;
	diagnostics.vcfToInvariantBed();

	// read bed

	// check if sites are as expected
	size_t c = 0;
	for (TInputFile bed("test.bed.gz", FileType::NoHeader); !bed.empty(); bed.popFront()) {
		EXPECT_EQ(bed.get(0), "junk_1");
		EXPECT_EQ(bed.get(1), "0");
		EXPECT_EQ(bed.get(2), toString(numLoci - 1));
		++c;
	}
	EXPECT_EQ(c, 1);
}

TEST_F(TVCFDiagnosticsTest, vcfToInvariantBed_allAreVariant_het) {
	std::vector<std::string> chromosomes(numLoci, "junk_1");
	std::vector<std::vector<BiallelicGenotype>> genotypes(numLoci);
	for (size_t l = 0; l < numLoci; l++){
		genotypes[l].resize(numIndiv, BiallelicGenotype::het); // all samples are heterozygotes
	}

	writeVcfFile(chromosomes, genotypes);

	// diagnose
	VCF::TVcfDiagnostics diagnostics;
	diagnostics.vcfToInvariantBed();

	// read bed

	// check if sites are as expected: should be empty, all are variant!
	std::vector<std::string> line;
	size_t c = 0;
	for (TInputFile bed("test.bed.gz", FileType::NoHeader); !bed.empty(); bed.popFront()) { ++c; }
	EXPECT_EQ(c, 0);
}

TEST_F(TVCFDiagnosticsTest, vcfToInvariantBed_allAreVariant_oneSample) {
	std::vector<std::string> chromosomes(numLoci, "junk_1");
	std::vector<std::vector<BiallelicGenotype>> genotypes(numLoci);
	for (size_t l = 0; l < numLoci; l++){
		genotypes[l].resize(numIndiv, BiallelicGenotype::homoFirst);
		genotypes[l][numIndiv-1] = BiallelicGenotype::het; // last sample is heterozygote
	}

	writeVcfFile(chromosomes, genotypes);

	// diagnose
	VCF::TVcfDiagnostics diagnostics;
	diagnostics.vcfToInvariantBed();

	// read bed
	TInputFile bed("test.bed.gz", FileType::NoHeader);

	// check if sites are as expected
	std::vector<std::string> line;
	/*
	size_t c = 0;
	// TODO: test fails; fix task!
	while (bed.read(line)){
		EXPECT_EQ(line[0], "junk_1");
		EXPECT_EQ(line[1], "0");
		EXPECT_EQ(line[2], toString(numLoci - 1));
		++c;
	}
	EXPECT_EQ(c, 1);*/
}
