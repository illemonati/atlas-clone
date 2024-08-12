#include "TGLFReader.h"
#include "genometools/VCF/TVcfWriter.h"

namespace GLF {

using genometools::Ploidy;

bool TGLFReader::_readChr() {
	// store position in file before read chromosome info
	auto posBeforeChr = _reader->tell() - 1; // substract record type

	// read chromosome info
	// name
	uint32_t len;
	if (!_reader->fill(len)) {
		_eof = true;
		return false;
	}

	std::string name;
	name.resize(len);
	_reader->fill(name);

	uint32_t refId;
	_reader->fill(refId);

	uint32_t length;
	_reader->fill(length);

	uint8_t ploidy;
	_reader->fill(ploidy);

	// check if chromosome is in index and set cur
	if(_hasIndex){
		_index.checkChromosome(refId, name, length, ploidy);
	} else {
		_index.addChromosme(name, length, ploidy, posBeforeChr);
	}

	// set first position = 0
	_position.move(refId, 0);
	_genotypeLikelihoodsGLF.type = _index.chromosomes()[refID()].isHaploid() ? Ploidy::haploid : Ploidy::diploid;

	return true;
};

bool TGLFReader::_readRecordType() {
	uint8_t tmpInt8;
	if (!_reader->fill(tmpInt8)) {
		_eof = true;
		return false;
	}
	_recordType = tmpInt8 >> 4;
	if (_recordType > 1) UERROR("Unknown record type in file '", _reader->name(), "'!");
	return true;
};

void TGLFReader::_readSNPRecord() {
	// read data of a single position
	// offset
	uint32_t offset;
	_reader->fill(offset);
	_position += offset;
	if (_position.position() >= _index.length(_position.refID())) {
		UERROR("Read a position in GLF file ", _reader->name(), " on chromosome ", refID(), " outside its length of ", _index.length(refID()), "!");
	}

	// maxLL and depth (both uint16)
	_reader->fill(_depth);
	_reader->fill(_RMS_mappingQual);

	// genotype likelihoods
	_reader->fill(_genotypeLikelihoodsGLF);
};


void TGLFReader::open(const std::string &Filename, bool HasIndex) {
	if (_reader->isOpen()) UERROR(_reader->name(), " is already open!");
	_reader.reset(coretools::makeReader(Filename));

	if (!_reader) UERROR("Failed to open GLF file '", Filename, "' for reading!");
	_index.clear();

	// parse header
	// version
	std::string fileVersion;
	fileVersion.resize(4);
	_reader->fill(fileVersion);

	if (fileVersion != version())
		UERROR("Non-supported GLF version '", fileVersion, "! Are you using GLF files produced with an earlier version of ATLAS?");

	// header
	uint32_t headerLen;
	_reader->fill(headerLen);
	if (headerLen > 0) {
		std::string header;
		header.resize(headerLen);
		if (_reader->fill(header) != headerLen) UERROR("Cannot read file ", Filename, "!");
	}
	_eof = false;

	// read index file
	_hasIndex = HasIndex;
	if(_hasIndex){
		_index.readChromosomes(Filename);
	}

	// read info of first chromosome
	popFront();
	//jumpToChr(0, false);
};

const genometools::TChromosomes& TGLFReader::chromosomes() const{
	assert(_hasIndex);
	return _index.chromosomes();
}

const genometools::TChromosome &TGLFReader::curChromosome() const {
	assert(_hasIndex);
	return _index.chromosomes()[refID()];
}

size_t TGLFReader::lastRefID() const {
	assert(_hasIndex);
	return _index.size() - 1;
}

void TGLFReader::rewind() {
	_eof = false;
	jumpToChr(0, false);
}

void TGLFReader::popFront() {
	// read record type
	if (!_readRecordType()) return;

	if (_recordType == 0) {
		_readChr();
		popFront();
	} else if (_recordType == 1) {
		_readSNPRecord();
	} else {
		UERROR("Unknown record type in file '", _reader->name(), "'!");
	}
};

bool TGLFReader::jumpToChr(uint32_t RefID, bool OnlyForward) {
	if (OnlyForward && refID() >= RefID) return true;
	_reader->seek(_index.positionInFile(RefID));

	// read record type
	if (!_readRecordType()){
		return false;
	}
	if(_recordType != 0){
		UERROR("File index (position in file) of chromosome '", _index.chromosomes()[RefID].name(), "' does not point to a valid position within the GLF file!");
	}

	_readChr();
	popFront();
	return !empty();
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
		popFront();
		if (empty()) return false;
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
		popFront();
		if (empty()) break;
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
	std::cout << curChr().name() << "\t" << _position.position() + 1 << "\t" << _depth << "\t" << (int)_RMS_mappingQual;
	for (size_t i = 0; i < _index.chrNumLikelihoodValues(refID()); ++i) std::cout << "\t" << _genotypeLikelihoodsGLF.data()[i];
	std::cout << "\n";
};

void TGLFReader::printToEnd() { // For debugging
	std::cout << "GLF version is " << version() << "\n";

	std::string oldChr = "";
	for(; !empty(); popFront()) {
		if (oldChr != curChr().name()) {
			printChr();
			oldChr = curChr().name();
		}
		printSite();
	}
}
void TGLFReader::writeIndex() {
	// read until end
	while (!empty()) {popFront();}

	//write index
	_index.writeChromosmes(_reader->name());
};
}
