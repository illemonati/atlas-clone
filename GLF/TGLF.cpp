/*
 * TGLF.cpp
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#include "TGLF.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <utility>
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "TGenotypeData.h"
#include "coretools/Types/probability.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Types/strongTypes.h"
#include "coretools/Types/weakTypes.h"

namespace GLF {
using namespace GenotypeLikelihoods;

//----------------------------------------------------
// TGlfChromosome
//----------------------------------------------------

namespace impl {

constexpr bool isHaploid(uint8_t Ploidy) {
	if (Ploidy == 1) {
		return true;
	} else if (Ploidy == 2) {
		return false;
	}
	UERROR("Currently GLFs only support ploidies 1 and 2 (not ", Ploidy, ")!");
};
} // namespace impl

TGlfChromosome::TGlfChromosome(const std::string &Name, uint32_t Length, uint8_t Ploidy)
    : _name(Name), _length(Length), _isHaploid(impl::isHaploid(Ploidy)) {};


void TGlfChromosome::update(std::string_view Name, uint16_t RefId, uint32_t Length, uint8_t Ploidy) {
	_name   = Name;
	_refId  = RefId;
	_length = Length;
	_isHaploid = impl::isHaploid(Ploidy);
};

void TGlfChromosome::clear() {
	_name      = "";
	_refId     = 0;
	_length    = 0;
	_isHaploid = false;
};

//---------------------------------
// TGlfWriter
//---------------------------------

void TGlfWriter::_writeHeader() {
	_write(_version.c_str(), 4 * sizeof(char));

	if (_header.length() > 0) {
		const uint32_t labelLength = _header.size();
		_write(&labelLength, sizeof(uint32_t));
		_write(_header.c_str(), labelLength * sizeof(char));
	} else {
		uint32_t zero32 = 0;
		_write(&zero32, sizeof(uint32_t));
	}
};

void TGlfWriter::open(const std::string &Filename, const std::string &Header) {
	_filename = Filename;
	_gzfp     = nullptr;
	_gzfp     = gzopen(_filename.c_str(), "wb");

	if (_gzfp == nullptr) UERROR("Failed to open file '", _filename, "' for writing!");
	_curChr.clear();

	// write header
	_header = Header;
	_writeHeader();
};

void TGlfWriter::newChromosome(const genometools::TChromosome &chromosome) {
	const uint8_t zero8 = 0;
	_write(&zero8, sizeof(uint8_t));

	// save cur info
	_curChr.update(chromosome.name(), chromosome.refID(), chromosome.length(), chromosome.ploidy());

	// write new chromosome: length of label, label, refId, length of ref sequence, ploidy
	const uint32_t labelLength = _curChr.name().size();
	_write(labelLength);
	_write(_curChr.name().c_str(), _curChr.name().length() * sizeof(char));
	_write(_curChr.refId());
	_write(_curChr.length());
	_write(_curChr.ploidy()); // TODO: I get an "uninitialized variable" error with valgrind. Why?

	// set oldPos and curChr
	_oldPos = 0;
};

void TGlfWriter::writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual,
			   GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods) {
	using genometools::Genotype;
	using coretools::Probability;
	const uint8_t _recordType1 = 1 << 4;
	// record type
	// TODO: add reference?
	_write(&_recordType1, sizeof(uint8_t));

	// offset
	const uint32_t offset = pos - _oldPos;
	_oldPos = pos;
	_write(&offset, sizeof(uint32_t));

	TGLFLikelihoods glfValues; // tmp used for writing

	// calculate likelihoods in GLF format
	// Note: genotype likelihoods are given for the 10 diploid genotypes!!
	if (_curChr.isHaploid()) {
		using genometools::Base;
		const double maxLik = std::max({genotypeLikelihoods[Genotype::AA], genotypeLikelihoods[Genotype::CC],
										genotypeLikelihoods[Genotype::GG], genotypeLikelihoods[Genotype::TT]});

		// normalize and scale to uint16
		glfValues.type = Ploidy::haploid;
		glfValues[Base::A] = Probability(genotypeLikelihoods[Genotype::AA] / maxLik);
		glfValues[Base::C] = Probability(genotypeLikelihoods[Genotype::CC] / maxLik);
		glfValues[Base::G] = Probability(genotypeLikelihoods[Genotype::GG] / maxLik);
		glfValues[Base::T] = Probability(genotypeLikelihoods[Genotype::TT] / maxLik);
	} else {
		// ploidy is 2
		glfValues.type = Ploidy::diploid;
		const double maxLik = *std::max_element(genotypeLikelihoods.begin(), genotypeLikelihoods.end());

		// normalize and scale to genometools::HighPrecisionPhredIntProbability

		for (auto g = Genotype::min; g < Genotype::max; ++g) {
			glfValues[g] = Probability(genotypeLikelihoods[g]/maxLik);
		}
	}

	// write maxLL as uint16_t
	// uint16_t maxLL_int = converter.toGlfFormat(maxLL);
	// write(&maxLL_int, sizeof(uint16_t));

	// write depth as uint16_t
	depth = std::min(depth, uint32_t{65535});
	const uint16_t tmp = depth;
	_write(&tmp, sizeof(uint16_t));

	// root mean square of mapping qualities
	_write(&RMS_mappingQual, sizeof(uint8_t));

	// genotype likelihoods
	_write(glfValues.data(), _curChr.numLikelihoodValues() * sizeof(genometools::HighPrecisionPhredIntProbability));
};

//---------------------------------
// TGlfReader
//---------------------------------
TGlfChromosome *TGlfReader::pointerToChr(uint32_t refId) {
	if (_curChr.refId() == refId) { return &_curChr; }
	
	auto it = _chromosomesAlreadyParsed.find(refId);
	return it == _chromosomesAlreadyParsed.end() ? nullptr : &it->second;
};

bool TGlfReader::_readChr() {
	// store current chromosome name in list of chromosomes parsed
	if (!_curChr.name().empty()) _chromosomesAlreadyParsed.emplace(_curChr.refId(), _curChr);

	// read chromosome info
	// name
	uint32_t len;
	if (!_read(&len, sizeof(uint32_t))) {
		_eof = true;
		return false;
	}
	constexpr size_t Nmax = 1024;
	char name[Nmax];
	_read(name, len);

	const auto refId  = _read<uint32_t>();
	const auto length = _read<uint32_t>();
	const auto ploidy = _read<uint8_t>();

	_curChr.update(std::string_view(name, len), refId, length, ploidy);

	// set first position = 0
	_position = 0;

	return true;
};

bool TGlfReader::_readRecordType() {
	uint8_t tmpInt8;
	if (!_read(&tmpInt8, sizeof(uint8_t))) {
		_eof = true;
		return false;
	}
	_recordType = tmpInt8 >> 4;
	if (_recordType > 1) UERROR("Unknown record type in file '", _filename, "'!");
	return true;
};

void TGlfReader::_readSNPRecord() {
	// read data of a single position
	// offset
	const auto offset = _read<uint32_t>();
	_position += offset;

	// maxLL and depth (both uint16)
	// read(&maxLL, sizeof(uint16_t));
	_depth           = _read<uint16_t>();
	_RMS_mappingQual = _read<uint8_t>();


	// genotype likelihoods
	_read(_genotypeLikelihoodsGLF.data(),
	      _curChr.numLikelihoodValues() * sizeof(genometools::HighPrecisionPhredIntProbability));
	_genotypeLikelihoodsGLF.type = _curChr.isHaploid() ? Ploidy::haploid : Ploidy::diploid;
};

void TGlfReader::setFilename(const std::string &Filename) { _filename = Filename; };

void TGlfReader::open(const std::string &Filename) {
	setFilename(Filename);
	_open();
};

void TGlfReader::_open() {
	if (_gzfp) UERROR(_filename, " is already open!");
	_gzfp = gzopen(_filename.c_str(), "rb");

	if (_gzfp == nullptr) UERROR("Failed to open file '", _filename, "' for reading!");
	_curChr.clear();
	_chromosomesAlreadyParsed.clear();
	_positionInFile = 0;

	// parse header
	// version
	char buffer[4];
	_read(buffer, 4 * sizeof(char));
	std::string fileVersion;
	fileVersion.assign(buffer, 4);

	if (fileVersion != _version)
		UERROR("Non-supported GLF version '", fileVersion, "! Currently only version '", _version, "' is supported.");
	_version = fileVersion;

	const auto headerLen = _read<uint32_t>();

	if (headerLen > 0) {
		char header[1024];
		if (_read(&header, headerLen) != (int) headerLen) UERROR("Cannot read file ", name(), "!");
	}
	_eof = false;

	// read info of first chromosome
	_readRecordType();
	if (_recordType != 0)
		UERROR("GLF file does not start with chromosome entry. The GLF format has changed with ATLAS 1.0, are you using "
			   "GLF files produced with an earlier version?");
	_readChr();
};

void TGlfReader::rewind() {
	// go back to beginning of file
	close();
	_open();
};

bool TGlfReader::readNext() {
	// read record type
	if (!_readRecordType()) return false;
	if (_recordType == 0) {
		_readChr();
		return readNext();
	} else if (_recordType == 1) {
		_readSNPRecord();
		return true;
	} else
		UERROR("Unknown record type in file '", _filename, "'!");
};

bool TGlfReader::_jumpToEndOfChr() {
	// read record type
	if (!_readRecordType()) return false;

	// tmp variables
	while (_recordType != 0) {
		// skipRecord();
		_readSNPRecord();

		if (!_readRecordType()) return false;
	}

	return true;
};

bool TGlfReader::jumpToNextChr() {
	if (!_jumpToEndOfChr()) { return false; }
	_readChr();
	return readNext();
};

bool TGlfReader::readNextWindow(std::vector<TGLFLikelihoods> &genoLikelihoods, uint32_t refId, uint32_t start,
				uint32_t end) {
	if (_eof) return false;

	// move to correct chromosome
	while (_curChr.refId() < refId) {
		if (!jumpToNextChr()) return false; // chromosome id doesn't exist (reached eof)
	}

	// jump to first position in window
	//  if there is a position 0 on the first chromosome parsed and windows starts right there
	//  -> we will not have read anything -> do this now
	while (_recordType == 0 || (_position < start && _curChr.refId() == refId)) {
		if (!readNext()) return false;
	}

	// have we passed window?
	if (_position >= end) return false; // no data

	// We are at first position in window with data

	// ensure size of container, fill with missing Data
	genoLikelihoods.resize(end - start, TGLFLikelihoods{});

	// Assumes that windows are read in order: no jumping back!
	if (refId < _curChr.refId()) { return false; }

	while (_position < end && _curChr.refId() == refId) {
		// fill in genotype likelihoods of current position
		genoLikelihoods[_position - start] = _genotypeLikelihoodsGLF;
		// read next record
		if (!readNext()) break;         // reached eof
	}
	return true;
};

// printing
void TGlfReader::printChr() {
	std::cout << "CHROMOSOME: '" << _curChr.name() << "' of length " << _curChr.length() << " and ploidy "
		  << (int)_curChr.ploidy() << "\n";
	
	// print header
	std::cout << "chr\tpos\tdepth\tRMS(MQ)";
	if(_curChr.isHaploid()){
		for(auto g = genometools::Base::min; g < genometools::Base::max; ++g){
			std::cout << "\tLL(" << genometools::toString(g) << ")";
		}		
	} else {
		for(auto g = genometools::Genotype::min; g < genometools::Genotype::max; ++g){
			std::cout << "\tLL(" << genometools::toString(g) << ")";
		}		
	}
	std::cout << "\n";
};

void TGlfReader::printSite() {
	// std::cout << curChr.name << "\t" << position << "\t" << maxLL << "\t" << depth << "\t" << RMS_mappingQual;
	//  print position as 1-based, internally it is 0-based
	std::cout << _curChr.name() << "\t" << _position + 1 << "\t" << _depth << "\t" << _RMS_mappingQual;
	for (size_t i = 0; i < _curChr.numLikelihoodValues(); ++i) std::cout << "\t" << _genotypeLikelihoodsGLF.data()[i];
	std::cout << "\n";
};

void TGlfReader::printToEnd() { // For debugging
	std::cout << "GLF version is " << _version << "\n";	

	// first chromosome	
	printChr();

	// now parse file
	std::string oldChr = _curChr.name();
	while (readNext()) {
		if (oldChr != _curChr.name()) {
			printChr();
			oldChr = _curChr.name();
		}
		printSite();
	}
};

}; // namespace GLF
