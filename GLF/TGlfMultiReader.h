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

#include "fmt/core.h"

#include "coretools/Containers/TBitSet.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Files/gzstream.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"

#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"

#include "TGLF.h"
#include "genometools/TFastaReader.h"

namespace GLF {

//----------------------------------------------------
// TMultiGLFDataSample
//----------------------------------------------------
class TMultiGLFDataSample {
private:
	TGLFLikelihoods _glf;//{{genometools::HighPrecisionPhredIntProbability::highest()}};
	size_t _depth = 0;
public:
	TMultiGLFDataSample(bool isHaploid = true)
		: _glf(genometools::HighPrecisionPhredIntProbability::highest(), Ploidy((!isHaploid))) {}
	constexpr TMultiGLFDataSample(const TGLFLikelihoods &GLs, uint16_t Depth) : _glf(GLs), _depth{Depth} {};

	constexpr bool hasData() const noexcept { return _depth > 0; };
	constexpr size_t depth() const noexcept { return _depth; };
	constexpr bool isHaploid() const noexcept { return _glf.isType(Ploidy::haploid); };

	constexpr genometools::HighPrecisionPhredIntProbability operator[](genometools::Genotype G) const noexcept {
		return _glf[G]; // asserts if diploid
	};

	constexpr genometools::HighPrecisionPhredIntProbability operator[](genometools::Base B) const noexcept {
		return _glf[B]; // asserts if diploid
	};
};

//-------------------------------------
// TGenotypeLikelihoodsOneAllelicCombination
//-------------------------------------
class TMultiGLFDataSampleOneAllelicCombination {
private:
	enum : uint8_t {NOTMISSING_DIPLOID = 0, MISSING_DIPLOID = 1, NOTMISSING_HAPLOID = 2, MISSING_HAPLOID = 3};

	coretools::TBitSet<2> _flags{MISSING_DIPLOID};
	std::array<coretools::Probability, 3> _GLs{coretools::Probability::highest(), coretools::Probability::highest(), coretools::Probability::highest()};

public:
	constexpr TMultiGLFDataSampleOneAllelicCombination(bool isHaploid=false) : _flags(MISSING_DIPLOID | isHaploid * MISSING_HAPLOID){}

	constexpr TMultiGLFDataSampleOneAllelicCombination(genometools::HighPrecisionPhredIntProbability homoFirst,
							   genometools::HighPrecisionPhredIntProbability het,
							   genometools::HighPrecisionPhredIntProbability homoSecond)
		: _flags(NOTMISSING_DIPLOID), _GLs({coretools::Probability(homoFirst), coretools::Probability(het), coretools::Probability(homoSecond)})  {}

	constexpr TMultiGLFDataSampleOneAllelicCombination(genometools::HighPrecisionPhredIntProbability first,
							   genometools::HighPrecisionPhredIntProbability second)
		:  _flags(NOTMISSING_HAPLOID), _GLs({coretools::Probability(first), coretools::Probability(second), coretools::Probability::highest()})  {}

	constexpr bool isMissing() const noexcept { return _flags.get<0>();}
	constexpr bool isHaploid() const noexcept { return _flags.get<1>();}

	constexpr coretools::Probability
	operator[](genometools::BiallelicGenotype Genotype) const noexcept {
		assert(isHaploid() == genometools::isHaploid(Genotype));
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
	coretools::TOutputFile _vcf;

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
			_vcf.writeNoDelim(genotypeStrings[mleGeno], ":0:");
		} else {
			auto slq = _getSecondHighestGTL(in, GTL);
			_vcf.writeNoDelim(genotypeStrings[mleGenotypes.front()], ':', genometools::PhredIntProbability(slq - minQual).get(), ':');
		}

		return minQual;
	}

	template<size_t NumGeno>
	void _writeCell(bool IsMissing, size_t Depth, const std::array<genometools::HighPrecisionPhredIntProbability, NumGeno> &GTL){
		if (IsMissing) {
			if constexpr (NumGeno == 3) { _vcf.write("./.:.:.:."); } // diploid
			else { _vcf.write(".:.:.:."); } // haploid
			return;
		}

		// write genotype and genotype quality
		const auto minQual = _writeGenotypeAndQuality(GTL);

		// write depth
		_vcf.writeNoDelim(Depth, ':');

		// write likelihoods
		_writeLikelihood(GTL[0] - minQual);
		for (size_t g = 1; g < NumGeno; g++){
			_vcf.writeNoDelim(',');
			_writeLikelihood(GTL[g] - minQual);
		}
		_vcf.writeDelim();
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
	bool _onlyPositionsWithData = false;
	std::vector<bool> _GLFIsActive;
	std::vector<TGlfReader *> _activeGLFs;

	// Moving along active files
	uint32_t _position = 0;
	uint32_t _curRefId = 0;
	TGlfChromosome _curChr;
	uint32_t _numActiveFilesWithData = 0;
	uint32_t _minDepth = 0;

	// reference
	genometools::TFastaReader fastaReader;

	void _openGLFs();
	int _getGLFIndexFromName(const std::string &name) const;
	void _setActive(size_t index);
	void _setAllInactive();
	void _prepareParsing();
	bool _jumpToNextPosition();

	bool _moveToNextChromosome();
	TMultiGLFData _data;

	static constexpr size_t _windowSize = 1000;
	size_t _iWindow = 0;
	std::vector<TMultiGLFData> _dataWindow;
	std::vector<size_t> _numActive;

public:
	//const TMultiGLFData& data() const noexcept {return _data;};
	//TMultiGLFData& data() noexcept {return _data;};
	const TMultiGLFData& data() const noexcept {return _dataWindow[_iWindow];};
	TMultiGLFData& data() noexcept {return _dataWindow[_iWindow];};

	TGlfMultiReader();
	TGlfMultiReader(const std::vector<std::string>& FileNames);

	~TGlfMultiReader();

	void openGLFs(const std::vector<std::string> &Filenames);
	void openGLFs();
	void closeGLF();
	void setDepthFilter(int MinDepth);
	void addReference(const std::string& FastaFile);
	void onlyPositionsWithData(bool set = true) { _onlyPositionsWithData = set; };

	// set active / inactive
	void setActive(int index);
	void setActive(const std::string &name);
	void setActive(int index1, int index2);
	void setActive(const std::string &name1, const std::string &name2);
	void setActive(std::vector<int> &indexes);
	void setActive(std::vector<std::string> &names);
	void setAllActive();

	// parse
	bool readWindow();
	bool readNextOld();
	bool readNext();

	// output
	std::vector<std::string> namesOfActiveFiles() const;
	std::vector<std::string> sampleNamesOfActiveFiles() const;

	// access data
	constexpr uint32_t numSamples() const noexcept { return _numGLFs; };
	uint32_t numActiveSamples() const noexcept { return _activeGLFs.size(); };
	//constexpr uint32_t numActiveSamplesWithData() const noexcept { return _numActiveFilesWithData; };
	constexpr uint32_t numActiveSamplesWithData() const noexcept { return _numActive[_iWindow]; };
	std::string chr() const { return _curChr.name(); };
	constexpr uint32_t position() const noexcept { return _position; };
	genometools::Base refBase() const noexcept {
		return fastaReader.isOpen() ? fastaReader(_curRefId, _position) : genometools::Base::N;
	};
};

}; // end namespace GLF

#endif /* GLF_TGLFMULTIREADER_H_ */
