/*
 * TGLF.cpp
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#include "TGLF.h"
#include "coretools/Main/TParameters.h"
#include "genometools/GenomePositions/TChromosomes.h"

namespace GLF {

bool TGLFReader::_readChr() {
	// store position in file before read chromosome info
	auto posBeforeChr = _positionInFile - 1; // substract record type

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

	// check if chromosome is in index and set cur
	if(_hasIndex){
		_index.checkChromosome(refId, std::string_view(name, len), length, ploidy);
	} else {
		_index.addChromosme(std::string_view(name, len), length, ploidy, posBeforeChr);
	}

	// set first position = 0
	_position.move(refId, 0);

	return true;
};

bool TGLFReader::_readRecordType() {
	uint8_t tmpInt8;
	if (!_read(&tmpInt8, sizeof(uint8_t))) {
		_eof = true;
		return false;
	}
	_recordType = tmpInt8 >> 4;
	if (_recordType > 1) UERROR("Unknown record type in file '", _filename, "'!");
	return true;
};

void TGLFReader::_readSNPRecord() {
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
		  _index.chrNumLikelihoodValues(refID()) * sizeof(genometools::HighPrecisionPhredIntProbability));

	_genotypeLikelihoodsGLF.type = _index.chromosomes()[refID()].isHaploid() ? Ploidy::haploid : Ploidy::diploid;
};

void TGLFReader::setFilename(const std::string &Filename) { _filename = Filename; };

void TGLFReader::open(const std::string &Filename, bool HasIndex) {
	setFilename(Filename);
	_open(HasIndex);
};

void TGLFReader::_open(bool HasIndex) {
	if (_gzfp) UERROR(_filename, " is already open!");
	_gzfp = gzopen(_filename.c_str(), "rb");

	if (_gzfp == nullptr) UERROR("Failed to open GLF file '", _filename, "' for reading!");
	_index.clear();
	_positionInFile = 0;

	// parse header
	// version
	char buffer[4];
	_read(buffer, 4 * sizeof(char));
	_glfFileVersion.assign(buffer, 4);

	if (_glfFileVersion != _version)
		UERROR("Non-supported GLF version '", _glfFileVersion, "! Are you using GLF files produced with an earlier version of ATLAS?");

	// header
	const auto headerLen = _read<uint32_t>();
	if (headerLen > 0) {
		char header[1024];
		if (_read(&header, headerLen) != (int) headerLen) UERROR("Cannot read file ", name(), "!");
	}
	_eof = false;

	// read index file
	_hasIndex = HasIndex;
	if(_hasIndex){
		_index.readChromosomes(_filename);
	}

	// read info of first chromosome
	_readRecordType();
	if (_recordType != 0)
		UERROR("GLF file does not start with chromosome entry. Are you using GLF files produced with an earlier version of ATLAS?");
	_readChr();
};

const genometools::TChromosomes& TGLFReader::chromosomes(){
	assert(_hasIndex);
	return _index.chromosomes();
}

const genometools::TChromosome &TGLFReader::curChromosome() {
	assert(_hasIndex);
	return _index.chromosomes()[refID()];
}

size_t TGLFReader::lastRefID() {
	assert(_hasIndex);
	return _index.size() - 1;
}

void TGLFReader::rewind() {
	// go back to beginning of file
	close();
	_open();
}

bool TGLFReader::readNext() {
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

bool TGLFReader::jumpToChr(uint32_t RefID, bool OnlyForward) {
	if (OnlyForward && refID() >= RefID) return true;
	if(gzseek64(_gzfp, _index.positionInFile(RefID), SEEK_SET) < 0){ return false; }

	// read record type
	if (!_readRecordType()){ return false; }
	if(_recordType != 0){
		UERROR("File index (position in file) of chromosome '", _index.chromosomes()[RefID].name(), "' does not point to a valid position within the GLF file!");
	}

	_readChr();
	return readNext();
};

bool TGLFReader::jumpToNextChr() {
	if(refID() == _index.size() - 1) { return false; }
	return jumpToChr(refID() + 1);
}
bool TGLFReader::jumpToPositionOrBeyond(const genometools::TGenomePosition &Position, bool OnlyForward) {
	// move to correct chromosome
	jumpToChr(Position.refID(), OnlyForward);


	// jump to first position at or after Position
	// Assume linear GLF access, i.e. if _position > Position, assume we are at correct position
	while (_recordType == 0 || _position < Position) {
		if (!readNext()) return false;
	}
	return true;
};

bool TGLFReader::readNextWindow(std::vector<TGLFLikelihoods> &GenoLikelihoods, genometools::TGenomeWindow Window) {
	if (_eof) return false;

	if(!jumpToPositionOrBeyond(Window.from())){ return false; }

	// have we passed window?
	if (_position > Window) return false; // no data

	// We are at first position in window with data

	// ensure size of container, fill with missing Data
	GenoLikelihoods.resize(Window.size(), TGLFLikelihoods{});

	// Assumes that windows are read in order: no jumping back!
	if (Window < _position) { return false; }

	while(Window.within(_position)){
		// fill in genotype likelihoods of current position
		GenoLikelihoods[_position - Window.from()] = _genotypeLikelihoodsGLF;
		// read next record
		if (!readNext()) break;         // reached eof
	}
	return true;
};

// printing
void TGLFReader::printChr() {
	std::cout << "CHROMOSOME: '" << curChr().name() << "' of length " << curChr().length() << " and ploidy "
			  << (int)curChr().ploidy() << "\n";

	// print header
	std::cout << "chr\tpos\tdepth\tRMS(MQ)";
	if(curChr().isHaploid()){
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

void TGLFReader::printSite() {
	// std::cout << curChr.name << "\t" << position << "\t" << maxLL << "\t" << depth << "\t" << RMS_mappingQual;
	//  print position as 1-based, internally it is 0-based
	std::cout << curChr().name() << "\t" << _position.position() + 1 << "\t" << _depth << "\t" << _RMS_mappingQual;
	for (size_t i = 0; i < _index.chrNumLikelihoodValues(refID()); ++i) std::cout << "\t" << _genotypeLikelihoodsGLF.data()[i];
	std::cout << "\n";
};

void TGLFReader::printToEnd() { // For debugging
	std::cout << "GLF version is " << _version << "\n";

	// first chromosome
	printChr();

	// now parse file
	std::string oldChr = curChr().name();
	while (readNext()) {
		if (oldChr != curChr().name()) {
			printChr();
			oldChr = curChr().name();
		}
		printSite();
	}
}
void TGLFReader::writeIndex() {
	// read until end
	while(readNext()){};

	//write index
	_index.writeChromosmes(_filename);
};

void TGLFPrinter::run() {
	TGLFReader reader(coretools::instances::parameters().get<std::string>("glf"));
	reader.printToEnd();
}

void TGLFIndexer::run(){
	// open GLF for reading, do not open index
	TGLFReader reader(coretools::instances::parameters().get<std::string>("glf"), false);
	reader.writeIndex();
}

}; // namespace GLF
