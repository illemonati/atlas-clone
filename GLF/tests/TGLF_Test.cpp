//
// Created by caduffm on 9/1/20.
//

#include "gtest/gtest.h"
#include <algorithm>

#include "genometools/GenotypeTypes.h"
#include "TGLF.h"
#include "TGenotypeData.h"
#include "TTestGLFFile.h"
#include "coretools/Strings/stringFunctions.h"

//-------------------------------------------------------------
// TGLF_Test_ReadWrite
// --> writing and reading
//-------------------------------------------------------------

using namespace GLF;
using genometools::Base;
using genometools::Genotype;
using genometools::HighPrecisionPhredIntProbability;
using GenotypeLikelihoods::TBaseLikelihoods;

class TGLF_Test_WriteRead : public ::testing::Test {
protected:
	std::string _filename = "testGLF.glf.gz";
	std::vector<size_t> chrLength;

public:
	GLF::TTestGLFFile outputGLF;
	GLF::TGlfReader inputGLF;

	// store stuff from input GLF to later compare with output
	std::vector<size_t> positions;
	std::vector<size_t> depths;
	std::vector<TGLFLikelihoods> genotypeLikelihoods;

	void write(int numDummySites) {
		// settings
		chrLength = {100, 200, 300};

		// open GLF file for writing
		outputGLF.openOutput(_filename, chrLength);

		// write sites
		outputGLF.writeDummySites(numDummySites);
		outputGLF.closeOutput();
	}

	void writeWithMissingSites() {
		// write some more complicated stuff: missing sites, depth = 0, mapping qual = 0, missing chromosomes etc.

		// settings
		chrLength = {50, 100, 150, 200, 250};

		// open GLF file for writing
		outputGLF.openOutput(_filename, chrLength);

		// write sites
		//  1) first chromosome is empty
		outputGLF.writeNewChromosome();

		// second chromosome has sites, but some sites have depth 0 or base N
		outputGLF.writeNewChromosome();
		// 2) ok site
		outputGLF.writeDummySite(10, 10);
		// 3) depth = 0
		outputGLF.writeDummySite(20, 0);
		// 4) depth = 10, but all bases are C
		std::vector<GenotypeLikelihoods::TBaseLikelihoods> bases;
		bases.reserve(10);
		Base base                    = Base::C;
		coretools::Probability error = 0.001;
		for (size_t d = 0; d < 10; d++) {
			TBaseLikelihoods baseData = GenotypeLikelihoods::fromError(base, error);
			bases.emplace_back(baseData);
		}
		GenotypeLikelihoods::TGenotypeLikelihoods gtL = GenotypeLikelihoods::getGLH(bases);
		outputGLF.writeDummySite(30, 0, gtL);
		// 5) depth = 10, all bases are A, but mapping quality is zero
		bases.clear();
		base = Base::A;
		for (size_t d = 0; d < 10; d++) {
			TBaseLikelihoods baseData = GenotypeLikelihoods::fromError(base, error);
			bases.emplace_back(baseData);
		}
		gtL = GenotypeLikelihoods::getGLH(bases);
		outputGLF.writeSite(40, 0, gtL, 0);

		// 6) third chromosome is empty
		outputGLF.writeNewChromosome();

		// 7) fourth chromosome has sites
		outputGLF.writeNewChromosome();
		// only last site of chromosome is not missing
		outputGLF.writeDummySite(199, 10);

		// 8) last chromosome is empty
		outputGLF.writeNewChromosome();

		outputGLF.closeOutput();
	}

	void writeWithDifferentPloidies() {
		// write chromosomes with different ploidies
		chrLength                     = {50, 100, 150};
		std::vector<uint8_t> ploidies = {2, 1, 2}; // second chromosome is haploid

		// open GLF file for writing
		outputGLF.openOutput(_filename, chrLength, ploidies);

		// write sites
		//  1) first chromosome is diploid
		outputGLF.writeNewChromosome();
		outputGLF.writeDummySite(10);
		outputGLF.writeDummySite(20);

		// 2) second chromosome is haploid
		outputGLF.writeNewChromosome();
		outputGLF.writeDummySite(10);
		outputGLF.writeDummySite(20);

		// 3) third chromosome is diploid
		outputGLF.writeNewChromosome();
		outputGLF.writeDummySite(10);
		outputGLF.writeDummySite(20);

		outputGLF.closeOutput();
	}
	virtual void read(bool open=true) {
		// open GLF for reading
		if (open) inputGLF.open(_filename);

		// read!
		while (inputGLF.readNext()) {
			// store
			positions.push_back(inputGLF.position());
			depths.push_back(inputGLF.depth());
			genotypeLikelihoods.push_back(inputGLF.genotypeLikelihoodsGLF());
		}
	}

	void TearDown() override{};
};

TEST_F(TGLF_Test_WriteRead, _name) {
	write(300);
	read();
	EXPECT_EQ(inputGLF.name(), _filename);
}

TEST_F(TGLF_Test_WriteRead, positions) {
	write(300);
	read();
	// check if written and read positions are equal
	int c = 0;
	for (auto writtenPosition = outputGLF.beginPositions(); writtenPosition != outputGLF.endPositions();
	     writtenPosition++, c++) {
		EXPECT_EQ(writtenPosition->position(), positions[c]);
	}
}

TEST_F(TGLF_Test_WriteRead, depth) {
	write(300);
	read();
	// check if written and read depth are equal
	int c = 0;
	for (auto writtenDepth = outputGLF.beginDepths(); writtenDepth != outputGLF.endDepths(); writtenDepth++, c++) {
		EXPECT_EQ(*writtenDepth, (int)depths[c]);
	}
}

TEST_F(TGLF_Test_WriteRead, chromosomes) {
	write(300);
	read();
	// check if written and read chromosomes are equal
	for (int i = 0; i < 3; i++) {
		const auto chr = inputGLF.pointerToChr(i);

		EXPECT_EQ(chr->refId(), i);
		EXPECT_EQ(chr->name(), "Chr" + coretools::str::toString(i + 1));
		EXPECT_EQ(chr->length(), chrLength[i]);
		EXPECT_EQ(chr->isHaploid(), false);
		EXPECT_EQ(chr->numLikelihoodValues(), 10); // diploid
	}
}

void normalizeByMax_Diploid(GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods) {
	coretools::normalize(genotypeLikelihoods, *std::max_element(genotypeLikelihoods.begin(), genotypeLikelihoods.end()));
}

void normalizeByMax_Haploid(GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods) {
	// find maxLL
	double maxLL = genotypeLikelihoods[Genotype::AA];
	if (genotypeLikelihoods[Genotype::CC] > maxLL) maxLL = genotypeLikelihoods[Genotype::CC];
	if (genotypeLikelihoods[Genotype::GG] > maxLL) maxLL = genotypeLikelihoods[Genotype::GG];
	if (genotypeLikelihoods[Genotype::TT] > maxLL) maxLL = genotypeLikelihoods[Genotype::TT];

	// normalize
	genotypeLikelihoods[Genotype::AA].scale(maxLL);
	genotypeLikelihoods[Genotype::CC].scale(maxLL);
	genotypeLikelihoods[Genotype::GG].scale(maxLL);
	genotypeLikelihoods[Genotype::TT].scale(maxLL);
}

TEST_F(TGLF_Test_WriteRead, genotypeLikelihoods) {
	write(300);
	read();
	// check if written and read genotype likelihoods are equal
	int c = 0;
	for (auto writtenGTL = outputGLF.beginGenotypeLikelihoods(); writtenGTL != outputGLF.endGenotypeLikelihoods();
	     writtenGTL++, c++) {
		// need to normalize the written likelihoods by maximal LL in order to compare
		normalizeByMax_Diploid(*writtenGTL);
		for (Genotype g = Genotype::min; g < Genotype::max; ++g) { // go over all 10 possible genotypes
			// compare in GLF format (quite a large imprecision when going from likelihood -> GLF likelihood ->
			// likelihood)
			EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[g]), genotypeLikelihoods[c][g]);
		}
	}
}

TEST_F(TGLF_Test_WriteRead, rewind) {
	write(300);
	read();
	positions.clear();
	inputGLF.rewind();
	read(false);

	// check if written and read positions are equal
	int c = 0;
	for (auto writtenPosition = outputGLF.beginPositions(); writtenPosition != outputGLF.endPositions();
	     writtenPosition++, c++) {
		EXPECT_EQ(writtenPosition->position(), positions[c]);
	}
}

TEST_F(TGLF_Test_WriteRead, jumpToNextChr) {
	write(300);
	read();
	inputGLF.rewind();

	// we start on chromosome 1
	EXPECT_EQ(0, inputGLF.refId());
	EXPECT_EQ(100, inputGLF.chrLength());

	// jump to first site on chromosome 2
	EXPECT_TRUE(inputGLF.jumpToNextChr());
	EXPECT_EQ(1, inputGLF.refId());
	EXPECT_EQ(200, inputGLF.chrLength());
	size_t whichSite = 0;
	for (auto pos = outputGLF.beginPositions(); pos != outputGLF.endPositions(); pos++, whichSite++) {
		if (pos->refID() == 1) break;
	}
	EXPECT_EQ(outputGLF.position(whichSite).position(), inputGLF.position());
	EXPECT_EQ(outputGLF.depth(whichSite), inputGLF.depth());

	// jump to next chromosome (=3) -> skip all other sites of chromosome 2
	EXPECT_TRUE(inputGLF.jumpToNextChr());
	EXPECT_EQ(2, inputGLF.refId());
	EXPECT_EQ(300, inputGLF.chrLength());
	whichSite = 0;
	for (auto pos = outputGLF.beginPositions(); pos != outputGLF.endPositions(); pos++, whichSite++) {
		if (pos->refID() == 2) break;
	}
	EXPECT_EQ(outputGLF.position(whichSite).position(), inputGLF.position());
	EXPECT_EQ(outputGLF.depth(whichSite), inputGLF.depth());

	// chromosome 3 was the last one -> jumping to the next one returns false
	EXPECT_FALSE(inputGLF.jumpToNextChr());
}

//-------------------------------------------------------------
// Check writing and reading with missing data
//-------------------------------------------------------------

TEST_F(TGLF_Test_WriteRead, positions_missingData) {
	writeWithMissingSites();
	read();
	// check if written and read positions are equal
	int c = 0;
	for (auto writtenPosition = outputGLF.beginPositions(); writtenPosition != outputGLF.endPositions();
	     writtenPosition++, c++) {
		EXPECT_EQ(writtenPosition->position(), positions[c]);
	}
}

TEST_F(TGLF_Test_WriteRead, depth_missingData) {
	writeWithMissingSites();
	read();
	// check if written and read depth are equal
	int c = 0;
	for (auto writtenDepth = outputGLF.beginDepths(); writtenDepth != outputGLF.endDepths(); writtenDepth++, c++) {
		EXPECT_EQ(*writtenDepth, (int)depths[c]);
	}
}

TEST_F(TGLF_Test_WriteRead, chromosomes_missingData) {
	writeWithMissingSites();
	read();
	// check if written and read chromosomes are equal
	// first chromosome
	for (int i = 0; i < 5; i++) {
		const auto chr = inputGLF.pointerToChr(i);

		EXPECT_EQ(chr->refId(), i);
		EXPECT_EQ(chr->name(), "Chr" + coretools::str::toString(i + 1));
		EXPECT_EQ(chr->length(), chrLength[i]);
		EXPECT_EQ(chr->isHaploid(), false);
		EXPECT_EQ(chr->numLikelihoodValues(), 10); // diploid
	}
}

TEST_F(TGLF_Test_WriteRead, genotypeLikelihoods_missingData) {
	writeWithMissingSites();
	read();
	// check if written and read genotype likelihoods are equal
	int c = 0;
	for (auto writtenGTL = outputGLF.beginGenotypeLikelihoods(); writtenGTL != outputGLF.endGenotypeLikelihoods();
	     writtenGTL++, c++) {
		// need to normalize the written likelihoods by maximal LL in order to compare
		normalizeByMax_Diploid(*writtenGTL);
		for (Genotype g = Genotype::min; g < Genotype::max; ++g) { // go over all 10 possible genotypes
			// compare in GLF format (quite a large imprecision when going from likelihood -> GLF likelihood ->
			// likelihood)
			EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[g]), genotypeLikelihoods[c][g]);
		}
	}
}

TEST_F(TGLF_Test_WriteRead, jumpToNextChr_missingData) {
	writeWithMissingSites();
	read();
	positions.clear();
	depths.clear();
	inputGLF.rewind();

	// we start on chromosome 1
	EXPECT_EQ(0, inputGLF.refId());
	EXPECT_EQ(50, inputGLF.chrLength());

	// jump to first site on chromosome 2
	inputGLF.jumpToNextChr();
	EXPECT_EQ(10, inputGLF.position());
	EXPECT_EQ(10, inputGLF.depth());
	EXPECT_EQ(1, inputGLF.refId());
	EXPECT_EQ(100, inputGLF.chrLength());

	// jump to next chromosome (=3), but this is empty, so jump to chromosome 4 and reads first site on chr4 -> skip all
	// other sites of chromosome 2
	inputGLF.jumpToNextChr();
	EXPECT_EQ(3, inputGLF.refId());
	EXPECT_EQ(200, inputGLF.chrLength());
	EXPECT_EQ(199, inputGLF.position());
	EXPECT_EQ(10, inputGLF.depth());

	// jump to chromosome 5 -> this is the last one -> returns false
	EXPECT_FALSE(inputGLF.jumpToNextChr());
	EXPECT_EQ(4, inputGLF.refId());
	EXPECT_EQ(250, inputGLF.chrLength());
}

//-------------------------------------------------------------
// Check writing and reading with chromosomes of different ploidies
//-------------------------------------------------------------

TEST_F(TGLF_Test_WriteRead, positions_withDifferentPloidies) {
	writeWithDifferentPloidies();
	read();
	// check if written and read positions are equal
	int c = 0;
	for (auto p = outputGLF.beginPositions(); p != outputGLF.endPositions(); p++, c++) {
		EXPECT_EQ(p->position(), positions[c]);
	}
}

TEST_F(TGLF_Test_WriteRead, depth_withDifferentPloidies) {
	writeWithDifferentPloidies();
	read();
	// check if written and read depth are equal
	int c = 0;
	for (auto writtenDepth = outputGLF.beginDepths(); writtenDepth != outputGLF.endDepths(); writtenDepth++, c++) {
		EXPECT_EQ(*writtenDepth, (int)depths[c]);
	}
}

TEST_F(TGLF_Test_WriteRead, chromosomes_withDifferentPloidies) {
	writeWithDifferentPloidies();
	read();
	// check if written and read chromosomes are equal
	// first chromosome
	for (int i = 0; i < 3; i++) {
		const auto chr = inputGLF.pointerToChr(i);

		EXPECT_EQ(chr->refId(), i);
		EXPECT_EQ(chr->name(), "Chr" + coretools::str::toString(i + 1));
		EXPECT_EQ(chr->length(), chrLength[i]);
		if (i == 1) { // haploid
			EXPECT_EQ(chr->isHaploid(), true);
			EXPECT_EQ(chr->numLikelihoodValues(), 4);
		} else { // diploid
			EXPECT_EQ(chr->isHaploid(), false);
			EXPECT_EQ(chr->numLikelihoodValues(), 10);
		}
	}
}

TEST_F(TGLF_Test_WriteRead, genotypeLikelihoods_withDifferentPloidies) {
	writeWithDifferentPloidies();
	read();
	// check if written and read genotype likelihoods are equal
	int c = 0;
	for (auto writtenGTL = outputGLF.beginGenotypeLikelihoods(); writtenGTL != outputGLF.endGenotypeLikelihoods();
	     writtenGTL++, c++) {
		if (c == 2 || c == 3) {
			normalizeByMax_Haploid(*writtenGTL);
			// haploid -> compare homozygous genotypes! ATTENTION: if haploid, only 4 values are returned -> access with
			// ACGT, not AA CC GG TT!!
			EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[Genotype::AA]), genotypeLikelihoods[c][Base::A]);
			EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[Genotype::GG]), genotypeLikelihoods[c][Base::G]);
			EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[Genotype::TT]), genotypeLikelihoods[c][Base::T]);
			EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[Genotype::CC]), genotypeLikelihoods[c][Base::C]);
		} else {
			// need to normalize the written likelihoods by maximal LL in order to compare
			normalizeByMax_Diploid(*writtenGTL);
			for (Genotype g = Genotype::min; g < Genotype::max; ++g) { // go over all 10 possible genotypes
				// compare in GLF format (quite a large imprecision when going from likelihood -> GLF likelihood ->
				// likelihood)
				EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[g]), genotypeLikelihoods[c][g]);
			}
		}
	}
}

//-------------------------------------------------------------
// TGLF_Test_WriteRead_Windows
// -- read in windows
//-------------------------------------------------------------

class TGLF_Test_WriteRead_Windows : public TGLF_Test_WriteRead {
protected:
	static constexpr size_t windowLen = 20;
	std::vector<std::vector<TGLFLikelihoods>> genotypeLikelihoods_perWindow;
public:
	void read(bool open=true) override {
		// open GLF for reading
		if (open) inputGLF.open(_filename);

		// parse GLFs in windows
		while (!inputGLF.eof()) {
			// move to new chromosome
			const auto curRefId    = inputGLF.refId();
			const auto curChrLen   = inputGLF.chrLength();
			size_t windowStart = 0;
			size_t windowEnd   = windowLen;

			// parse all windows of chromosome (unless file ends before chromosome ends)
			while (windowStart < curChrLen && !inputGLF.eof()) {
				// resize storage
				std::vector<TGLFLikelihoods> genoLikelihoods_oneWindow;

				// read data
				const auto isGood = inputGLF.readNextWindow(genoLikelihoods_oneWindow, curRefId, windowStart, windowEnd);
				if (isGood || inputGLF.eof()) {
					// store if window contains data
					genotypeLikelihoods_perWindow.push_back(genoLikelihoods_oneWindow);
				}
				// move window
				windowStart = windowEnd;
				windowEnd   = windowStart + windowLen;
			}
		}
	}

	void TearDown() override {}
};

TEST_F(TGLF_Test_WriteRead_Windows, genotypeLikelihoods_writeAll) {
	write(600);
	read();
	// check if written and read genotype likelihoods are equal
	auto writtenGTL = outputGLF.beginGenotypeLikelihoodsWithMissingSites();
	for (auto &window : genotypeLikelihoods_perWindow) {
		for (auto genotypeLikelihood_read : window) {
			// need to normalize the written likelihoods by maximal LL in order to compare
			normalizeByMax_Diploid(*writtenGTL);
			for (Genotype g = Genotype::min; g < Genotype::max; ++g) { // go over all 10 possible genotypes
				// compare in GLF format (quite a large imprecision when going from likelihood -> GLF likelihood ->
				// likelihood)
				EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[g]), genotypeLikelihood_read[g]);
			}
			writtenGTL++;
		}
	}
}

TEST_F(TGLF_Test_WriteRead_Windows, genotypeLikelihoods_writeWithMissingSites) {
	write(300);
	read();
	// check if written and read genotype likelihoods are equal
	auto writtenGTL = outputGLF.beginGenotypeLikelihoodsWithMissingSites();
	for (auto &window : genotypeLikelihoods_perWindow) {
		for (auto genotypeLikelihood_read : window) {
			// need to normalize the written likelihoods by maximal LL in order to compare
			normalizeByMax_Diploid(*writtenGTL);
			if (genotypeLikelihood_read.type == Ploidy::diploid) {
				for (Genotype g = Genotype::min; g < Genotype::max; ++g) { // go over all 10 possible genotypes
					// compare in GLF format (quite a large imprecision when going from likelihood -> GLF likelihood ->
					// likelihood)
					EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[g]), genotypeLikelihood_read[g]);
				}
			} else {
				for (Base b = Base::min; b < Base::max; ++b) { // go over all 4 possible genotypes
					EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[genotype(b, b)]),
							  genotypeLikelihood_read[b]);
				}
			}
			writtenGTL++;
		}
	}
}

TEST_F(TGLF_Test_WriteRead_Windows, genotypeLikelihoods_writeWithMissingWindows) {
	writeWithMissingSites();
	read();

	EXPECT_EQ(genotypeLikelihoods_perWindow.size(), 4); // only 4 windows with data

	// first window: chromosome 2, 0-19
	auto writtenGTL = outputGLF.beginGenotypeLikelihoodsWithMissingSites() + 50;
	for (int s = 0; s < 20; s++) {
		normalizeByMax_Diploid(*writtenGTL);
		if (genotypeLikelihoods_perWindow[0][s].type == Ploidy::diploid) {
			for (Genotype g = Genotype::min; g < Genotype::max; ++g) { // go over all 10 possible genotypes
				EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[g]), genotypeLikelihoods_perWindow[0][s][g]);
			}
		} else {
			for (Base b = Base::min; b < Base::max; ++b) { // go over all 10 possible genotypes
				EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[genotype(b, b)]),
						  genotypeLikelihoods_perWindow[0][s][b]);
			}
		}
		writtenGTL++;
	}

	// second window: chromosome 2, 20-39
	writtenGTL = outputGLF.beginGenotypeLikelihoodsWithMissingSites() + 70;
	for (int s = 0; s < 20; s++) {
		normalizeByMax_Diploid(*writtenGTL);
		if (genotypeLikelihoods_perWindow[1][s].type == Ploidy::diploid) {
			for (Genotype g = Genotype::min; g < Genotype::max; ++g) { // go over all 10 possible genotypes
				EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[g]), genotypeLikelihoods_perWindow[1][s][g]);
			}
		} else {
			for (Base b = Base::min; b < Base::max; ++b) { // go over all 10 possible genotypes
				EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[genotype(b, b)]),
						  genotypeLikelihoods_perWindow[1][s][b]);
			}
		}

		writtenGTL++;
	}

	// third window: chromosome 2, 40-59
	writtenGTL = outputGLF.beginGenotypeLikelihoodsWithMissingSites() + 90;
	for (int s = 0; s < 20; s++) {
		normalizeByMax_Diploid(*writtenGTL);
		if (genotypeLikelihoods_perWindow[2][s].type == Ploidy::diploid) {
			for (Genotype g = Genotype::min; g < Genotype::max; ++g) { // go over all 10 possible genotypes
				EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[g]), genotypeLikelihoods_perWindow[2][s][g]);
			}
		} else {
			for (Base b = Base::min; b < Base::max; ++b) { // go over all 10 possible genotypes
				EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[genotype(b, b)]),
						  genotypeLikelihoods_perWindow[2][s][b]);
			}
		}
		writtenGTL++;
	}

	// fourth window: chromosome 4, 180-199
	writtenGTL = outputGLF.beginGenotypeLikelihoodsWithMissingSites() + 480;
	for (int s = 0; s < 20; s++) {
		normalizeByMax_Diploid(*writtenGTL);
		if (genotypeLikelihoods_perWindow[3][s].type == Ploidy::diploid) {
			for (Genotype g = Genotype::min; g < Genotype::max; ++g) { // go over all 10 possible genotypes
				EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[g]), genotypeLikelihoods_perWindow[3][s][g]);
			}
		} else {
			for (Base b = Base::min; b < Base::max; ++b) { // go over all 10 possible genotypes
				EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[genotype(b, b)]),
						  genotypeLikelihoods_perWindow[3][s][b]);
			}
		}
		writtenGTL++;
	}
}

TEST_F(TGLF_Test_WriteRead_Windows, oneWindow_writeAll) {
	write(200);

	// open GLF for reading
	inputGLF.open(_filename);

	// resize storage
	std::vector<TGLFLikelihoods> genoLikelihoods_oneWindow;
	inputGLF.readNextWindow(genoLikelihoods_oneWindow, 1, 100, 120);
	genotypeLikelihoods_perWindow.push_back(genoLikelihoods_oneWindow);

	// check if written and read genotype likelihoods are equal
	auto writtenGTL = outputGLF.beginGenotypeLikelihoodsWithMissingSites() + 200;
	for (auto &window : genotypeLikelihoods_perWindow) {
		for (auto genotypeLikelihood_read : window) {
			// need to normalize the written likelihoods by maximal LL in order to compare
			normalizeByMax_Diploid(*writtenGTL);
			if (genotypeLikelihood_read.type == Ploidy::diploid) {
				for (Genotype g = Genotype::min; g < Genotype::max; ++g) { // go over all 10 possible genotypes
					// compare in GLF format (quite a large imprecision when going from likelihood -> GLF likelihood ->
					// likelihood)
					EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[g]), genotypeLikelihood_read[g]);
				}
			} else {
				for (Base b = Base::min; b < Base::max; ++b) { // go over all 10 possible genotypes
					// compare in GLF format (quite a large imprecision when going from likelihood -> GLF likelihood ->
					// likelihood)
					EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[genotype(b, b)]),
							  genotypeLikelihood_read[b]);
				}
			}
			writtenGTL++;
		}
	}
}

TEST_F(TGLF_Test_WriteRead_Windows, oneWindow_writeWithMissingSites) {
	writeWithMissingSites();

	// open GLF for reading
	inputGLF.open(_filename);

	// resize storage
	std::vector<TGLFLikelihoods> genoLikelihoods_oneWindow(windowLen);
	inputGLF.readNextWindow(genoLikelihoods_oneWindow, 3, 180, 200);
	genotypeLikelihoods_perWindow.push_back(genoLikelihoods_oneWindow);

	// check if written and read genotype likelihoods are equal
	// fourth window: chromosome 4, 180-199
	auto writtenGTL = outputGLF.beginGenotypeLikelihoodsWithMissingSites() + 480;
	for (int s = 0; s < 20; s++) {
		normalizeByMax_Diploid(*writtenGTL);
		if (genotypeLikelihoods_perWindow[0][s].type == Ploidy::diploid) {
			for (Genotype g = Genotype::min; g < Genotype::max; ++g) { // go over all 10 possible genotypes
				EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[g]), genotypeLikelihoods_perWindow[0][s][g]);
			}
		} else {
			for (Base b = Base::min; b < Base::max; ++b) { // go over all 10 possible genotypes
				EXPECT_EQ(HighPrecisionPhredIntProbability((*writtenGTL)[genotype(b, b)]), genotypeLikelihoods_perWindow[0][s][b]);
			}
		}
		writtenGTL++;
	}
}
