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

#include "coretools/Containers/TStrongArray.h"
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

	TMultiGLFDataSampleOneAllelicCombination(genometools::HighPrecisionPhredIntProbability homoFirst,
							   genometools::HighPrecisionPhredIntProbability het,
							   genometools::HighPrecisionPhredIntProbability homoSecond)
		: _flags(NOTMISSING_DIPLOID), _GLs({coretools::Probability(homoFirst), coretools::Probability(het), coretools::Probability(homoSecond)})  {}

	TMultiGLFDataSampleOneAllelicCombination(genometools::HighPrecisionPhredIntProbability first,
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

TMultiGLFDataOneAllelicCombination fill(coretools::TConstView<GLF::TMultiGLFDataSample> samples,
										genometools::AllelicCombination alleleicCombination);

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
	size_t _minSamplesWithData  = 0;
	std::vector<bool> _GLFIsActive;
	std::vector<TGlfReader *> _activeGLFs;

	// Moving along active files
	TGlfChromosome _curChr;
	uint32_t _curRefId = 0;
	uint32_t _minDepth = 0;
	size_t _windowStart = 0;
	size_t _windowSize  = 64;
	std::vector<TMultiGLFData> _dataWindow;
	std::vector<size_t> _numActive;

	// reference
	genometools::TFastaReader fastaReader;

	void _openGLFs();
	int _getGLFIndexFromName(const std::string &name) const;
	void _setActive(size_t index);
	void _setAllInactive();
	void _prepareParsing();
	bool _jumpToNextPosition();

	bool _moveToNextChromosome();

public:
	const TMultiGLFData& data(size_t iWindow) const noexcept {return _dataWindow[iWindow];};

	TGlfMultiReader();

	void openGLFs(const std::vector<std::string> &Filenames);
	void openGLFs();
	void addReference(const std::string& FastaFile);
	void minSamplesWithData(size_t MinSamplesWithData) { _minSamplesWithData = MinSamplesWithData; };

	// set active / inactive
	void setActive(int index);
	void setActive(const std::string &name);
	void setActive(int index1, int index2);
	void setActive(const std::string &name1, const std::string &name2);
	void setActive(const std::vector<int> &indexes);
	void setActive(const std::vector<std::string> &names);
	void setAllActive();

	// parse
	std::vector<size_t> readWindow();

	// output
	std::vector<std::string> namesOfActiveFiles() const;
	std::vector<std::string> sampleNamesOfActiveFiles() const;

	// access data
	uint32_t numActiveSamples() const noexcept { return _activeGLFs.size(); }
	std::string chr() const { return _curChr.name(); }
	constexpr uint32_t position(size_t iWindow) const noexcept { return _windowStart + iWindow; }
	bool hasRef() const noexcept {return fastaReader.isOpen();}
	genometools::Base refBase(size_t iWindow) const noexcept {
		return hasRef() ? fastaReader(_curRefId, position(iWindow)) : genometools::Base::N;
	}
};

}; // end namespace GLF

#endif /* GLF_TGLFMULTIREADER_H_ */
