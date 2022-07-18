/*
 * TGlfMultiReader.h
 *
 *  Created on: Feb 11, 2020
 *      Author: wegmannd
 */

#ifndef GLF_TGLFMULTIREADER_H_
#define GLF_TGLFMULTIREADER_H_

#include <stddef.h>
#include <stdint.h>
#include <algorithm>
#include <array>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TBitSet.h"
#include "TFastaBuffer.h"
#include "TGLF.h"
#include "TGenomePosition.h"
#include "TRandomGenerator.h"
#include "gzstream.h"

namespace GLF {

//----------------------------------------------------
// TMultiGLFDataSample
//----------------------------------------------------
class TMultiGLFDataSample {
private:
	const TGLFLikelihoods *_glf;
	union {
		uint16_t _depth;
		bool _isHaploid;
	};
public:
	TMultiGLFDataSample() = default;
	constexpr TMultiGLFDataSample(bool isHaploid) : _glf(nullptr), _isHaploid{isHaploid} {}
	constexpr TMultiGLFDataSample(const TGLFLikelihoods *GLs, uint16_t Depth) : _glf(GLs), _depth{Depth}{};

	constexpr bool hasData() const noexcept { return _glf; };
	constexpr uint16_t depth() const noexcept { return hasData() * _depth; };
	constexpr bool isHaploid() const noexcept { return hasData() ? _glf->type == Ploidy::haploid : _isHaploid; };

	template<bool HasData>
	constexpr bool isHaploid() const noexcept {
		if constexpr (HasData) return _glf->type == Ploidy::haploid;
		else return _isHaploid;
	}

	constexpr const genometools::HighPrecisionPhredIntProbability get_with_data(const genometools::Genotype &G) const noexcept {
		return (*_glf)[G]; // asserts if haploid
	};

	constexpr const genometools::HighPrecisionPhredIntProbability get_with_data(const genometools::Base &B) const noexcept {
		return (*_glf)[B];// asserts if diploid
	};

	constexpr const genometools::HighPrecisionPhredIntProbability operator[](const genometools::Genotype &G) const noexcept {
		if (!_glf) { return genometools::HighPrecisionPhredIntProbability::highest(); }
		return get_with_data(G);
	};

	constexpr const genometools::HighPrecisionPhredIntProbability operator[](const genometools::Base &B) const noexcept {
		if (!_glf) { return genometools::HighPrecisionPhredIntProbability::highest(); }
		return get_with_data(B);
	};
};

//-------------------------------------
// TGenotypeLikelihoodsOneAllelicCombination
//-------------------------------------
class TMultiGLFDataSampleOneAllelicCombination {
private:
	enum : uint8_t {NOTMISSING_DIPLOID = 0, MISSING_DIPLOID = 1, NOTMISSING_HAPLOID = 2, MISSING_HAPLOID = 3};

	coretools::TBitSet<2> _flags{MISSING_DIPLOID};
	std::array<genometools::HighPrecisionPhredIntProbability, 3> _GLs{genometools::HighPrecisionPhredIntProbability::highest(), genometools::HighPrecisionPhredIntProbability::highest(), genometools::HighPrecisionPhredIntProbability::highest()};

public:
	constexpr TMultiGLFDataSampleOneAllelicCombination(bool isHaploid=false) : _flags(MISSING_DIPLOID | isHaploid * MISSING_HAPLOID){}

	constexpr TMultiGLFDataSampleOneAllelicCombination(genometools::HighPrecisionPhredIntProbability homoFirst,
							   genometools::HighPrecisionPhredIntProbability het,
							   genometools::HighPrecisionPhredIntProbability homoSecond)
		: _flags(NOTMISSING_DIPLOID), _GLs({homoFirst, het, homoSecond})  {}

	constexpr TMultiGLFDataSampleOneAllelicCombination(genometools::HighPrecisionPhredIntProbability first,
							   genometools::HighPrecisionPhredIntProbability second)
		:  _flags(NOTMISSING_HAPLOID), _GLs({first, second, genometools::HighPrecisionPhredIntProbability::highest()})  {}

	constexpr bool isMissing() const noexcept { return _flags.get<0>();}
	constexpr bool isHaploid() const noexcept { return _flags.get<1>();}

	constexpr genometools::HighPrecisionPhredIntProbability
	operator[](genometools::BiallelicGenotype Genotype) const noexcept {
		//assert(isHaploid() && !genometools::isHaploid(Genotype));
		//assert(!isHaploid() && genometools::isHaploid(Genotype));
		//if (isMissing()) { return genometools::HighPrecisionPhredIntProbability::highest(); }
		return _GLs[genometools::altAlleleCounts(Genotype)];
	};
};

using TMultiGLFDataOneAllelicCombination = std::vector<TMultiGLFDataSampleOneAllelicCombination>;
using TMultiGLFData                      = std::vector<TMultiGLFDataSample>;

void fill(TMultiGLFDataOneAllelicCombination &storage, const TMultiGLFData &samples,
	  genometools::AllelicCombination alleleicCombination);

//----------------------------------------------------
// TGlfMultiReaderVcfLocusDefinition
//----------------------------------------------------
class TGlfMultiReaderVcf {
private:
	bool _usePhredScaledLikelihoods;

protected:
	mutable gz::ogzstream _vcf;

	genometools::Base _ref{};
	genometools::Base _alt{};

	genometools::Genotype _refHom{};
	genometools::Genotype _het{};
	genometools::Genotype _altHom{};

	void _openVCF(const std::string & Filename, const std::string & Source, const std::vector<std::string> & SampleNames);
	void _closeVCF();
	void _setMajorMinor(genometools::Base refAllele, genometools::Base altAllele);
	void _writeLikelihood(genometools::HighPrecisionPhredIntProbability likGlf);
	void _writeSiteInformation(const std::string & ChrName, uint32_t Position, genometools::PhredIntProbability VariantQuality, size_t Depth);

	template<size_t NumGeno>
	auto _getSecondHighestGTL(size_t IndexBest,
	                          const std::array<genometools::HighPrecisionPhredIntProbability, NumGeno> &GTL) {
		if constexpr (NumGeno == 2) { // haploid
			return GTL[1 - IndexBest]; // the other one
		} else {
			constexpr std::array<std::array<size_t, 2>, NumGeno> a = {{{1, 2}, {0, 2}, {0, 1}}}; // the two other ones
			return std::min(GTL[a[IndexBest][0]], GTL[a[IndexBest][1]]);
		}
	}

	template<size_t NumGeno>
	auto _writeGenotypeAndQuality(const std::array<genometools::HighPrecisionPhredIntProbability, NumGeno> &GTL) {
		std::array<std::string, NumGeno> genotypeStrings;
		if constexpr (NumGeno == 2){ genotypeStrings = {"0", "1"}; } // haploid
		else { genotypeStrings = {"0/0", "0/1", "1/1"}; } // diploid

		const auto min     = std::min_element(GTL.cbegin(), GTL.cend());
		const auto minQual = *min;
		const auto in      = std::distance(GTL.cbegin(), min);

		// get all genotypes with minQual (=MLE)
		std::vector<size_t> mleGenotypes;
		for (size_t i = 0; i < NumGeno; ++i) {
			if (GTL[i] == minQual) { mleGenotypes.push_back(i); }
		}

		// write genotype quality
		if (mleGenotypes.size() > 1) {
			const auto mleGeno = mleGenotypes[coretools::instances::randomGenerator().sample(mleGenotypes.size())];
			_vcf << '\t' << genotypeStrings[mleGeno] << ':';
			_vcf << "0:"; // difference between best and second best is 0
		} else {
			_vcf << '\t' << genotypeStrings[mleGenotypes.front()] << ':';
			// find second highest quality
			auto slq = _getSecondHighestGTL(in, GTL);
			_vcf << genometools::PhredIntProbability(slq - minQual) << ":";
		}

		return minQual;
	}

	template<size_t NumGeno>
	void _writeCell(bool IsMissing, size_t Depth, const std::array<genometools::HighPrecisionPhredIntProbability, NumGeno> &GTL){
		if (IsMissing) {
			if constexpr (NumGeno == 3) { _vcf << "\t./.:.:.:."; } // diploid
			else { _vcf << "\t.:.:.:."; } // haploid
			return;
		}

		// write genotype and genotype quality
		const auto minQual = _writeGenotypeAndQuality(GTL);

		// write depth
		_vcf << Depth << ':';

		// write likelihoods
		_writeLikelihood(GTL[0] - minQual);
		for (size_t g = 1; g < NumGeno; g++){
			_vcf << ',';
			_writeLikelihood(GTL[g] - minQual);
		}
	}

public:
	TGlfMultiReaderVcf(const std::string & Filename, const std::string & Source, const std::vector<std::string> &sampleNames,
			   bool UsePhredScaledLikelihoods);
	~TGlfMultiReaderVcf() { _closeVCF(); }

	void writeSite(const std::string &chrName, uint32_t position, genometools::PhredIntProbability varianTQuality,
		       TMultiGLFData &data, genometools::Base Ref, genometools::Base Alt);
};

//----------------------------------------------------
// TGlfMultiReader
//----------------------------------------------------
class TGlfMultiReader {
private:
	size_t _numGLFs = 0;
	std::vector<std::string> _GLFNames;
	std::vector<TGlfReader> _GLFs;
	bool _readersOpened = false;

	// active files
	// Object will loop only over active files
	bool _onlyJumpToPositionsWithData = false;
	std::vector<bool> _GLFIsActive;
	std::vector<TGlfReader *> _pointerToActiveGLFs;

	// Moving along active files
	uint32_t _position = 0;
	uint32_t _nextPosition = 0; // next is anticipated position, used to advance
	uint32_t _curRefId = 0;
	TGlfChromosome* _curChr;
	uint32_t _numActiveFilesWithData = 0;
	uint32_t _minDepth = 0;

	// reference
	bool hasReference = false;
	BAM::TFastaBuffer fastaBuffer;

	void _openGLFs();
	int _getGLFIndexFromName(const std::string &name) const;
	void _setActive(size_t index);
	void _setAllInactive();
	void _prepareParsing();
	bool _jumpToNextPosition();

	bool _moveToNextChromosome();

	void _writeDiploidIndividualToVCF(int ind, gz::ogzstream &vcf, genometools::Base major, genometools::Base minor,
					 const std::vector<std::string> &genotypeStrings, bool usePhredLikelihoods);
	void _writeHaploidIndividualToVCF(int ind, gz::ogzstream &vcf, genometools::Base major, genometools::Base minor,
					 const std::vector<std::string> &genotypeStrings, bool usePhredLikelihoods);

public:
	TMultiGLFData data;

	TGlfMultiReader();
	TGlfMultiReader(const std::vector<std::string>& FileNames);

	~TGlfMultiReader();

	void openGLFs(const std::vector<std::string> &Filenames);
	void openGLFs();
	void closeGLF();
	void setDepthFilter(int MinDepth);
	void addReference(const std::string& FastaFile);
	void onlyJumpToPositionsWithData(bool set = true) { _onlyJumpToPositionsWithData = set; };

	// set active / inactive
	void setActive(int index);
	void setActive(const std::string &name);
	void setActive(int index1, int index2);
	void setActive(const std::string &name1, const std::string &name2);
	void setActive(std::vector<int> &indexes);
	void setActive(std::vector<std::string> &names);
	void setAllActive();

	// parse
	bool readNext();

	// output
	void print() const;
	void writeSampleNamesOfActiveFiles(gz::ogzstream &out, std::string& sep) const;
	std::vector<std::string> namesOfActiveFiles() const;
	std::vector<std::string> sampleNamesOfActiveFiles() const;

	// access data
	constexpr uint32_t numSamples() const noexcept { return _numGLFs; };
	uint32_t numActiveSamples() const noexcept { return _pointerToActiveGLFs.size(); };
	constexpr uint32_t numActiveSamplesWithData() const noexcept { return _numActiveFilesWithData; };
	std::string chr() const { return _curChr->name(); };
	constexpr uint32_t position() const noexcept { return _position; };
	constexpr genometools::Base refBase() const noexcept {
		return hasReference ? fastaBuffer.refAt(genometools::TGenomePosition(_curRefId, _position)) : genometools::Base::N;
	};
};

}; // end namespace GLF

#endif /* GLF_TGLFMULTIREADER_H_ */
