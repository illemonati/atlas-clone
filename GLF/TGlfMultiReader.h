/*
 * TGlfMultiReader.h
 *
 *  Created on: Feb 11, 2020
 *      Author: wegmannd
 */

#ifndef GLF_TGLFMULTIREADER_H_
#define GLF_TGLFMULTIREADER_H_

#include "GenotypeTypes.h"
#include "TBitSet.h"
#include "TFastaBuffer.h"
#include "TGLF.h"
#include "stringFunctions.h"

namespace GLF {

//----------------------------------------------------
// TMultiGLFDataSample
//----------------------------------------------------

class TMultiGLFDataSample {
private:
	const TGLFLikelihoods *_glf = nullptr; // points to data TGlfReader
	uint16_t _depth_or_haploid  = 0;
public:
	constexpr TMultiGLFDataSample() = default;
	constexpr TMultiGLFDataSample(bool isHaploid) : _depth_or_haploid(isHaploid) {}
	constexpr TMultiGLFDataSample(const TGLFLikelihoods *GLs, uint16_t Depth) : _glf(GLs), _depth_or_haploid(Depth){};

	constexpr bool hasData() const noexcept { return _glf != nullptr; };
	constexpr uint16_t depth() const noexcept { return hasData() ? _depth_or_haploid : 0; };
	constexpr bool isHaploid() const noexcept { return hasData() ? _glf->isHaploid() : _depth_or_haploid; };

	constexpr const genometools::HighPrecisionPhredIntProbability operator[](const genometools::Genotype &G) const {
		if (!hasData()) { return genometools::HighPrecisionPhredIntProbability::highest(); }
		return (*_glf)[G]; // throws if haploid
	};

	constexpr const genometools::HighPrecisionPhredIntProbability operator[](const genometools::Base &B) const {
		if (!hasData()) { return genometools::HighPrecisionPhredIntProbability::highest(); }
		return (*_glf)[B];// throws if diploid
	};
};

//-------------------------------------
// TGenotypeLikelihoodsOneAllelicCombination
//-------------------------------------
class TMultiGLFDataSampleOneAllelicCombination {
private:
	bool _isMissing = true;
	bool _isHaploid = false;
	std::array<genometools::HighPrecisionPhredIntProbability, 3> _GLs;

public:
	constexpr TMultiGLFDataSampleOneAllelicCombination(bool isHaploid=false) : _isHaploid(isHaploid) {}

	constexpr TMultiGLFDataSampleOneAllelicCombination(genometools::HighPrecisionPhredIntProbability homoFirst,
							   genometools::HighPrecisionPhredIntProbability het,
							   genometools::HighPrecisionPhredIntProbability homoSecond)
	    : _isMissing(false), _isHaploid(false), _GLs({homoFirst, het, homoSecond}) {}

	constexpr TMultiGLFDataSampleOneAllelicCombination(genometools::HighPrecisionPhredIntProbability first,
							   genometools::HighPrecisionPhredIntProbability second)
	    : _isMissing(false), _isHaploid(true), _GLs({first, second, genometools::HighPrecisionPhredIntProbability{}}) {}

	bool isMissing() const { return _isMissing; };
	bool isHaploid() const { return _isHaploid; };

	constexpr genometools::HighPrecisionPhredIntProbability
	operator[](genometools::BiallelicGenotype Genotype) const {
		// check
		if (_isHaploid && !genometools::isHaploid(Genotype)) {
		    throw std::runtime_error((std::string) __PRETTY_FUNCTION__ + ": sample is haploid!");
		}
		if (!_isHaploid && !genometools::isDiploid(Genotype)) {
		    throw std::runtime_error((std::string) __PRETTY_FUNCTION__ + ": sample is diploid!");
		}

		// return
		if (_isMissing) { return genometools::HighPrecisionPhredIntProbability::highest(); }
		return _GLs[genometools::altAlleleCounts(Genotype)];
	};
};

using TMultiGLFDataOneAllelicCombination = std::vector<TMultiGLFDataSampleOneAllelicCombination>;
using TMultiGLFData                      = std::vector<TMultiGLFDataSample>;

void fill(TMultiGLFDataOneAllelicCombination &storage, const TMultiGLFData &samples,
	  const genometools::AllelicCombination &alleleicCombination);

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
		if constexpr (NumGeno == 2) {
			return GTL[1 - IndexBest]; // the other one
		} else {
			constexpr std::array<std::array<size_t, 2>, NumGeno> a = {{{1, 2}, {0, 2}, {0, 1}}}; // the two other ones
			return std::min(GTL[a[IndexBest][0]], GTL[a[IndexBest][1]]);
		}
	}

	template<size_t NumGeno>
	auto _writeGenotypeAndQuality(const std::array<genometools::HighPrecisionPhredIntProbability, NumGeno> &GTL) {
		std::array<std::string, NumGeno> genotypeStrings;
		if constexpr (NumGeno == 2){ genotypeStrings = {"0", "1"}; }
		else { genotypeStrings = {"0/0", "0/1", "1/1"}; }

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
			if constexpr (NumGeno == 2) { _vcf << "\t./.:.:.:."; }
			else { _vcf << "\t.:.:.:."; }
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
	uint8_t numGLFs = 0;
	std::vector<std::string> GLFNames;
	std::vector<TGlfReader> GLFs;
	bool readersOpened = false;

	// active files
	// Object will loop only over active files
	bool _onlyJumpToPositionsWithData = false;
	std::vector<bool> GLFIsActive;
	std::vector<TGlfReader *> pointerToActiveGLFs;

	// Moving along active files
	uint32_t _position = 0;
	uint32_t _nextPosition = 0; // next is anticipated position, used to advance
	uint32_t _curRefId = 0;
	TGlfChromosome _curChr;
	uint32_t _numActiveFilesWithData = 0;
	uint32_t minDepth = 0;

	// reference
	bool hasReference = false;
	BAM::TFastaBuffer fastaBuffer;

	void _openGLFs();
	int _getGLFIndexFromName(const std::string &name) const;
	void _setActive(const int index);
	void _setAllInactive();
	void _prepareParsing();
	bool _jumpToNextPositionWithData();
	void _updateChromosomeInfo();


	bool moveToNextChromosome();

	void writeDiploidIndividualToVCF(int ind, gz::ogzstream &vcf, genometools::Base major, genometools::Base minor,
					 const std::vector<std::string> &genotypeStrings, bool usePhredLikelihoods);
	void writeHaploidIndividualToVCF(int ind, gz::ogzstream &vcf, genometools::Base major, genometools::Base minor,
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
	constexpr uint32_t numSamples() const noexcept { return numGLFs; };
	uint32_t numActiveSamples() const noexcept { return pointerToActiveGLFs.size(); };
	constexpr uint32_t numActiveSamplesWithData() const noexcept { return _numActiveFilesWithData; };
	std::string chr() const { return _curChr.name(); };
	constexpr uint32_t position() const noexcept { return _position; };
	constexpr genometools::Base refBase() const noexcept {
		return hasReference ? fastaBuffer.refAt(genometools::TGenomePosition(_curRefId, _position)) : genometools::Base::N;
	};
};

}; // end namespace GLF

#endif /* GLF_TGLFMULTIREADER_H_ */
