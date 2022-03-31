/*
 * TGlfMultireader.cpp
 *
 *  Created on: Feb 11, 2020
 *      Author: wegmannd
 */

#include "TGlfMultiReader.h"
#include "GenotypeTypes.h"
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "debugtools.h"
#include <algorithm>
#include <iterator>
#include <numeric>

namespace GLF {
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::instances::randomGenerator;
using genometools::HighPrecisionPhredIntProbability;

using coretools::str::toString;

void fill(TMultiGLFDataOneAllelicCombination &storage, const TMultiGLFData &samples,
	  const genometools::AllelicCombination &alleleicCombination) {
	using namespace genometools;
	storage.clear();
	storage.reserve(samples.size());
	for (const auto &s : samples) {
		if (s.isHaploid()) {
			if (s.hasData()) {
				storage.emplace_back(s[first(alleleicCombination)], s[second(alleleicCombination)]);
			} else {
				storage.emplace_back(true);
			}
		} else {
			if (s.hasData()) {
				storage.emplace_back(s[homoFirst(alleleicCombination)], s[het(alleleicCombination)],
						     s[homoSecond(alleleicCombination)]);
			} else {
				storage.emplace_back(false);
			}
		}
	}
};

uint32_t totalDepth(const TMultiGLFData &samples) noexcept {
	return std::accumulate(samples.begin(), samples.end(), uint32_t{}, [](auto tot, auto s) {return tot + s.depth();});
};

//----------------------------------------------------
// TGlfMultiReaderVcf
//----------------------------------------------------
TGlfMultiReaderVcf::TGlfMultiReaderVcf(const std::string filename, const std::string source,
				       std::vector<std::string> &sampleNames, bool UsePhredScaledLikelihoods)
    : _usePhredScaledLikelihoods(UsePhredScaledLikelihoods) {
	_openVCF(filename, source, sampleNames);
};

void TGlfMultiReaderVcf::_openVCF(const std::string &filename, const std::string &source,
				  std::vector<std::string> &sampleNames) {
	_closeVCF();

	// open vcf file
	vcf.open(filename.c_str());

	// write info
	// TODO: create VCF class to harmonize code across different uses. Also include code in Tiger and other
	vcf << "##fileformat=VCFv4.2\n";
	vcf << "##source=" << source << "\n";

	// make sure the header matches the format used in writeSiteToVCF
	vcf << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
	vcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	vcf << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n";
	vcf << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
	if (_usePhredScaledLikelihoods) {
		vcf << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled normalized genotype likelihoods\">\n";
	} else {
		vcf << "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Normalized genotype likelihoods\">\n";
	}

	// also write header with sample names
	vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	for (std::string &name : sampleNames) vcf << '\t' << name;
	vcf << '\n';
};

void TGlfMultiReaderVcf::_closeVCF() {
	vcf.close();
};

void TGlfMultiReaderVcf::_setMajorMinor(genometools::Base refAllele, genometools::Base altAllele) {
	_ref    = refAllele;
	_alt    = altAllele;
	_refHom = genometools::genotype(_ref, _ref);
	_het    = genometools::genotype(_ref, _alt);
	_altHom = genometools::genotype(_alt, _alt);
};

void TGlfMultiReaderVcf::writeLikelihood(HighPrecisionPhredIntProbability likGlf) const {
	if (_usePhredScaledLikelihoods) {
		vcf << (genometools::PhredIntProbability)likGlf;
	} else {
		if (likGlf == 0)
			vcf << '0';
		else
			vcf << (coretools::Log10Probability)likGlf;
	}
};

void TGlfMultiReaderVcf::writeSite(const std::string &chrName, uint32_t position,
				   const genometools::PhredIntProbability varianTQuality, TMultiGLFData &data,
				   genometools::Base refAllele, genometools::Base altAllele) {
	// Note: we pass hom/het indexes to maintain the major / minor order! Passing the alleleic combination is not enough
	// TODO: find way to harmonize code with TCaller
	_setMajorMinor(refAllele, altAllele);

	// write position
	vcf << chrName << '\t' << position + 1 << "\t.\t"; // internal position is zero-based!

	// write ref and alt alleles
	vcf << refAllele << '\t' << altAllele;

	// write quality of variant
	vcf << '\t' << varianTQuality;

	// write info field: total depth
	vcf << "\t.\tDP=" << totalDepth(data);

	// write filter, info and format
	if (_usePhredScaledLikelihoods)
		vcf << "\tGT:GQ:DP:PL";
	else
		vcf << "\tGT:GQ:DP:GL";

	// now write active samples
	for (const auto &d : data) {
		if (d.isHaploid())
			writeHaploidIndividualToVCF(d);
		else
			writeDiploidIndividualToVCF(d);
	}

	// end of line
	vcf << '\n';
};

void TGlfMultiReaderVcf::writeDiploidIndividualToVCF(const TMultiGLFDataSample &sample) const {
	if (!sample.hasData()) {
		vcf << "\t./.:.:.:.";
		return;
	}

	constexpr std::array genotypeStrings = {"0/0", "0/1", "1/1"};
	// find min qual
	const std::array smp = {sample[_refHom], sample[_het], sample[_altHom]};
	const auto min       = std::min_element(smp.cbegin(), smp.cend());
	const auto minQual   = *min;
	const auto in        = std::distance(smp.cbegin(), min);

	// get all genotypes with minQual (=MLE)
	std::vector<size_t> mleGenotypes;
	for (size_t i = 0; i < smp.size(); ++i)
		if (smp[i] == minQual) mleGenotypes.push_back(i);

	// write genotype quality
	if (mleGenotypes.size() > 1) {
		const int mleGeno = mleGenotypes[randomGenerator().sample(mleGenotypes.size())];
		vcf << '\t' << genotypeStrings[mleGeno] << ':';
		vcf << "0:";
	} else {
		vcf << '\t' << genotypeStrings[mleGenotypes.front()] << ':';
		constexpr std::array<std::array<size_t, 2>, 3> a = {{{1, 2}, {0, 2}, {1, 2}}};
		const auto slq = std::min(smp[a[in][0]], smp[a[in][1]]);
		// find second highest quality
		vcf << genometools::PhredIntProbability(slq - minQual) << ":";
	}

	// write depth
	vcf << sample.depth() << ':';

	// write likelihoods
	writeLikelihood(sample[_refHom] - minQual);
	vcf << ',';
	writeLikelihood(sample[_het] - minQual);
	vcf << ',';
	writeLikelihood(sample[_altHom] - minQual);
};

void TGlfMultiReaderVcf::writeHaploidIndividualToVCF(const TMultiGLFDataSample &sample) const {
	if (!sample.hasData()) {
		vcf << "\t./.:.:.:.";
		return;
	}

	constexpr std::array genotypeStrings = {"0", "1", "0/0", "0/1", "1/1"};
	using genometools::index;
	// find min qual
	const auto minQual = std::min(sample[_ref], sample[_alt]);

	// get all genotypes with minQual (=MLE)
	std::vector<genometools::Base> mleGenotypes;
	if (sample[_ref] == minQual) mleGenotypes.push_back(_ref);
	if (sample[_alt] == minQual) mleGenotypes.push_back(_alt);

	// write MLE genoytpe
	const genometools::Base mleGeno = mleGenotypes[randomGenerator().sample(mleGenotypes.size())];
	vcf << '\t' << genotypeStrings[index(mleGeno)] << ':';

	// write genotype quality
	if (mleGeno == _ref)
		vcf << (genometools::PhredIntProbability)(sample[_alt] - minQual) << ":";
	else
		vcf << (genometools::PhredIntProbability)(sample[_ref] - minQual) << ":";

	// write depth
	vcf << sample.depth() << ':';

	// write likelihoods
	writeLikelihood(sample[_ref] - minQual);
	vcf << ',';
	writeLikelihood(sample[_alt] - minQual);
};

//----------------------------------------------------
// TGlfMultiReader
//----------------------------------------------------
TGlfMultiReader::TGlfMultiReader() {
	setDepthFilter(parameters().getParameterWithDefault<int>("minDepth", 0));
};

	TGlfMultiReader::TGlfMultiReader(const std::vector<std::string>& FileNames) : GLFNames(FileNames){
	_openGLFs();
};

TGlfMultiReader::~TGlfMultiReader() {
	if (readersOpened) {
		closeGLF();
	}
};

void TGlfMultiReader::_openGLFs() {
	numGLFs     = GLFNames.size();
	GLFIsActive.resize(numGLFs, false);

	// open files
	GLFs.clear();
	GLFs.reserve(numGLFs);
	readersOpened = true;
	logfile().startIndent("Opening " + toString(numGLFs) + " GLF files:");
	for (const auto &name : GLFNames) {
		logfile().listFlush("Opening GLF '" + name + "' ...");
		GLFs.emplace_back(name);
		logfile().done();
	}
	logfile().endIndent();
	_setAllInactive();
};

void TGlfMultiReader::openGLFs(const std::vector<std::string> &FileNames) {
	GLFNames = FileNames;
	_openGLFs();
};

void TGlfMultiReader::openGLFs() {
	using namespace coretools::str;
	const auto parameter = parameters().getParameter<std::string>("glf");
	// assume that GLF file names are given in a file if string does not contain ".gz"
	if (!stringContains(parameter, ".gz")) {
		logfile().list("Reading glf input names from file '" + parameter + "'");
		std::ifstream in(parameter.c_str());
		if (!in) throw "Failed to open file '" + parameter + "'!";

		// read file
		while (in.good() && !in.eof()) {
			std::string line;
			std::getline(in, line);
			std::vector<std::string> vec;
			fillContainerFromStringWhiteSpace(line, vec, true);
			// skip empty lines
			if (vec.size() > 0) { GLFNames.push_back(vec[0]); }
		}
		in.close();
	} else {
		parameters().fillParameterIntoContainer("glf", GLFNames, ',');
	}
	_openGLFs();
};

void TGlfMultiReader::closeGLF() {
	if (readersOpened) {
		for (auto& g: GLFs) g.close();

		GLFNames.clear();
		GLFs.clear();
		pointerToActiveGLFs.clear();
		numGLFs       = 0;
		readersOpened = false;
	}
};

void TGlfMultiReader::setDepthFilter(int MinDepth) {
	minDepth = MinDepth;
	if (minDepth > 0) logfile().list("Will only keep sites with depth >= " + toString(minDepth) + ".");
};

void TGlfMultiReader::addReference(const std::string &FastaFile) {
	fastaBuffer.initialize(FastaFile, 1000000);
	hasReference = true;
};

//-------------------------------------
// set active / inactive
//-------------------------------------
int TGlfMultiReader::_getGLFIndexFromName(const std::string &name) const {
	const auto res = std::find(GLFNames.cbegin(), GLFNames.cend(), name);
	if (res == GLFNames.cend()) throw "GLF with name '" + name + "' not in TGlfMultiReader!";
	return std::distance(GLFNames.begin(), res);
};

void TGlfMultiReader::_setActive(const int index) {
	if (index >= numGLFs) throw "Index out of range in TGlfMultiReader::setActive(const int index)!";
	if (!GLFIsActive[index]) {
		GLFIsActive[index] = true;
		pointerToActiveGLFs.push_back(&GLFs[index]);
	}
};

void TGlfMultiReader::_setAllInactive() {
	GLFIsActive.assign(numGLFs, false);
	pointerToActiveGLFs.clear();
};

void TGlfMultiReader::_prepareParsing() {
	// prepare data
	data.clear();
	data.reserve(numActiveSamples());

	for (TGlfReader *it : pointerToActiveGLFs) it->rewind();

	// read first SNP record in all active files
	for (TGlfReader *it : pointerToActiveGLFs) { it->readNext(); }

	// where to start?
	if (_onlyJumpToPositionsWithData) {
		_jumpToNextPositionWithData();
	} else {
		_curRefId     = 0;
		_position     = 0;
		_nextPosition = 0;

		// fill chromosome info
		_updateChromosomeInfo();
	}
};

bool TGlfMultiReader::_jumpToNextPositionWithData() {
	// find min(refId) and min(pos)
	bool allFilesReachedEnd = true;
	bool first              = true;
	_curRefId               = pointerToActiveGLFs[0]->refId();
	_position               = pointerToActiveGLFs[0]->position();

	for (TGlfReader *it : pointerToActiveGLFs) {
		if (!it->eof()) {
			allFilesReachedEnd = false;

			if (first) {
				_curRefId = it->refId();
				_position = it->position();
				first     = false;
			} else {
				if (it->refId() < _curRefId) {
					_curRefId = it->refId();
					_position = it->position();
				} else if (it->refId() == _curRefId && it->position() < _position) {
					_position = it->position();
				}
			}
		}
	}

	// did we reach end?
	if (allFilesReachedEnd) return false;

	// fill chromosome info
	_updateChromosomeInfo();

	// make sure first record is used
	_nextPosition = _position;
	return true;
};

void TGlfMultiReader::_updateChromosomeInfo() {
	// update curChr
	TGlfChromosome *chr = pointerToActiveGLFs[0]->pointerToChr(_curRefId);
	if (chr == nullptr) {
		throw "Chromosome with refId " + toString(_curRefId) + " not present in GLF file '" +
		    pointerToActiveGLFs[0]->name() + "!";
	}
	_curChr = *chr;

	// check that all files share the same chromosome info
	chr   = nullptr;
	for (TGlfReader *it : pointerToActiveGLFs) {
		if (it->fillPointerToChr(_curRefId, chr)) {
			if (chr->name() != _curChr.name())
				throw "Chrosomome names differ between files '" + pointerToActiveGLFs[0]->name() + "' and '" +
				    it->name() + "': '" + _curChr.name() + "' != '" + chr->name() + "'!";
			if (chr->length() != _curChr.length())
				throw "Chrosomome lengths differ between files '" + pointerToActiveGLFs[0]->name() + "' and '" +
				    it->name() + "': '" + toString(_curChr.length()) + "' != '" + toString(chr->length()) + "'!";
		} else {
			throw "Chromosome with refId " + toString(_curRefId) +
			    " not present in any GLF file provided! Limit to only sites with data?";
		}
	}
};

void TGlfMultiReader::setActive(const int index) {
	_setAllInactive();
	_setActive(index);
	_prepareParsing();
};

void TGlfMultiReader::setActive(const std::string &name) {
	setActive(_getGLFIndexFromName(name));
};

void TGlfMultiReader::setActive(const int index1, const int index2) {
	_setAllInactive();
	_setActive(index1);
	_setActive(index2);
	_prepareParsing();
};

void TGlfMultiReader::setActive(const std::string &name1, const std::string &name2) {
	setActive(_getGLFIndexFromName(name1), _getGLFIndexFromName(name2));
};

void TGlfMultiReader::setActive(std::vector<int> &indexes) {
	_setAllInactive();
	for (const auto i : indexes) _setActive(i);
	_prepareParsing();
};

void TGlfMultiReader::setActive(std::vector<std::string> &names) {
	_setAllInactive();
	for (const auto &name : names) {
		_setActive(_getGLFIndexFromName(name));
	}
	_prepareParsing();
};

void TGlfMultiReader::setAllActive() {
	_setAllInactive();
	for (int i = 0; i < numGLFs; ++i) _setActive(i);
	_prepareParsing();
};

//-------------------------------------
// Looping over active files
//-------------------------------------
bool TGlfMultiReader::moveToNextChromosome() {
	// increment chromosome ref_char id
	++_curRefId;

	// advance all active files behind
	bool allFilesReachedEnd = true;
	for (TGlfReader *it : pointerToActiveGLFs) {
		while (!it->eof() && it->refId() < _curRefId) it->jumpToNextChr();
		if (!it->eof()) allFilesReachedEnd = false;
	}

	// check if we reached end of all files
	if (allFilesReachedEnd) return false;

	// get name and length from first active file not at end
	if (_onlyJumpToPositionsWithData) {
		_jumpToNextPositionWithData();
	} else {
		_position     = 0;
		_nextPosition = 0;
	}

	// get pointer to chromosome
	_updateChromosomeInfo();
	return true;
};

bool TGlfMultiReader::readNext() {
	// advance to next position
	if (_nextPosition > _curChr.length() && !moveToNextChromosome()) return false;

	// advance all files behind next position
	_numActiveFilesWithData = 0;
	bool allFilesReachedEnd = true;
	data.clear();
	for (TGlfReader *reader : pointerToActiveGLFs) {
		while (!reader->eof() && reader->refId() == _curRefId && reader->position() < _nextPosition) { reader->readNext(); }

		if (!reader->eof() && reader->position() == _nextPosition && reader->refId() == _curRefId) {
			data.emplace_back(&reader->genotypeLikelihoodsGLF(), reader->depth());
			++_numActiveFilesWithData;
		} else {
			data.emplace_back(reader->pointerToChr(_curRefId)->isHaploid()); // data is missing
		}
		if (!reader->eof()) allFilesReachedEnd = false;
	}

	// check if we reached end of all files
	if (allFilesReachedEnd) return false;

	// jump?
	if (_onlyJumpToPositionsWithData && _numActiveFilesWithData == 0) {
		if (_jumpToNextPositionWithData()) return readNext();
		return false;
	}

	// update position
	_position = _nextPosition;
	++_nextPosition;

	return true;
};

void TGlfMultiReader::print() const {
	std::cout << "\nMulti Reader at position " << _position << " on chromosome '" << _curChr.name() << std::endl;
	for (size_t i = 0; i < numActiveSamples(); ++i) {
		std::cout << "File " << i << ":";
		if (data[i].isHaploid()) {
			for (genometools::Base g = genometools::Base::min; g < genometools::Base::max; ++g) {
				std::cout << "\t" << data[i][g];
			}
		} else {
			for (genometools::Genotype g = genometools::Genotype::min; g < genometools::Genotype::max; ++g) {
				std::cout << "\t" << data[i][g];
			}
		}
		std::cout << std::endl;
	}
};

void TGlfMultiReader::writeSampleNamesOfActiveFiles(gz::ogzstream &out, std::string &sep) const {
	// sample names are file names without glf ending
	for (TGlfReader *it : pointerToActiveGLFs) { out << sep << coretools::str::readBeforeLast(it->name(), ".glf"); }
};

std::vector<std::string> TGlfMultiReader::namesOfActiveFiles() const {
	std::vector<std::string> vec;

	// sample names are file names without glf ending
	for (TGlfReader *it : pointerToActiveGLFs) { vec.emplace_back(it->name()); }
	return vec;
};

std::vector<std::string> TGlfMultiReader::sampleNamesOfActiveFiles() const {
	std::vector<std::string> vec;

	// sample names are file names without glf ending
	for (TGlfReader *it : pointerToActiveGLFs) { vec.emplace_back(coretools::str::readBeforeLast(it->name(), ".glf")); }
	return vec;
};

}; // end namespace GLF
