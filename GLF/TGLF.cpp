/*
 * TGLF.cpp
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#include "TGLF.h"
#include "coretools/Main/TParameters.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Files/TInputFile.h"


namespace GLF {
using namespace GenotypeLikelihoods;

//--------------------------------------------
// TGlfIndex
//--------------------------------------------

std::string TGlfIndex::_getIndexFileName(std::string_view FileName){
	std::string name(coretools::str::readBeforeLast(FileName, '.'));
	return name + ".idx";
}

void TGlfIndex::clear(){
	_chrs.clear();
	_posInFile.clear();
}

void TGlfIndex::addChromosme(std::string_view Name, uint32_t Length, uint8_t Ploidy, uint64_t PosInFile){
	if(Ploidy != 1 && Ploidy != 2)
		UERROR("Currently GLFs only support ploidies 1 and 2 (not ", Ploidy, ")!");
	_chrs.appendChromosome(Name, Length, Ploidy);
	_posInFile.push_back(PosInFile);
	_curChr = _posInFile.size() - 1;
}

void TGlfIndex::addChromosme(const genometools::TChromosome& Chr, uint64_t PosInFile){
	if(Chr.ploidy() != 1 && Chr.ploidy() != 2)
		UERROR("Currently GLFs only support ploidies 1 and 2 (not ", Chr.ploidy(), ")!");
	_chrs.appendChromosome(Chr);
	_posInFile.push_back(PosInFile);
	_curChr = _posInFile.size() - 1;
}

void TGlfIndex::readChromosomes(std::string_view GLFFilename) {
	clear();
	for(auto in = coretools::TInputFile(_getIndexFileName(GLFFilename), coretools::FileType::Header); in.isOpen(); in.popFront()){ 
		_chrs.appendChromosome(in.get("Chr"), in.get<uint32_t>("Length"), in.get<uint8_t>("Ploidy"));
		_posInFile.push_back(in.get<uint64_t>("PosInFile"));
	}
}

void TGlfIndex::writeChromosmes(std::string_view GLFFilename){
	coretools::TOutputFile out(_getIndexFileName(GLFFilename), {"Chr", "Length", "Ploidy", "PosInFile"});
	for(auto c: _chrs){
		out.writeln(c.name(), c.length(), c.ploidy(), _posInFile[c.refID()]);
	}
}

void TGlfIndex::jumpToNextChromosome() {
	++_curChr;
	if(_curChr == _chrs.size()){
		DEVERROR("Can not jump to next chromosome: already at end of GLF index!");
	}
}

void TGlfIndex::seekChr(uint32_t RefID) {
	if(RefID < _chrs.size()){
		_curChr = RefID;
	} else {
		UERROR("Chromosome with ref ID '", RefID, "' not found in GLF index!");
	}
}

void TGlfIndex::setCurChromosomeToAndCheck(std::string_view Name, uint32_t Length, uint8_t Ploidy, uint64_t PositionInFile) {
	if(!_chrs.exists(Name))
		UERROR("Chromosome '", Name, "' not found in GLF index!");
	_curChr = _chrs.refID(Name);

	if(_chrs[_curChr].length() != Length)
		UERROR("Length of chromosome '", Name, "' does not match length given in GLF index!");

	if(_chrs[_curChr].ploidy() != Ploidy)
		UERROR("Ploidy of chromosome '", Name, "' does not match ploidy given in GLF index!");

	if(_posInFile[_curChr] != PositionInFile)
		UERROR("File index (position in file) of chromosome '", Name, "' does not match value given in GLF index!");
}

bool TGlfIndex::hasSameChromosomes(const TGlfIndex &Other) const { 
	// checks if two TGlfIndexFiles contain the same chromosomes in terms of order, names and lengths
	// ploidy may be different (e.g. X for males and females), as well as posInFile, which depends on the amount of data
	if(_chrs.size() != Other._chrs.size()) return false;

	auto a = _chrs.begin();
	auto b = Other._chrs.begin();
	while(a != _chrs.end()){
		if(a->name() != b->name() || a->length() != b->length()) return false;
		++a; ++b;
	}
	return true;
}

//---------------------------------
// TGlfWriter
//---------------------------------

void TGlfWriter::_writeHeader(const std::string &Header) {
	_write(_version.c_str(), 4 * sizeof(char));
	
	_header = Header;
	if (_header.length() > 0) {
		const uint32_t labelLength = _header.size();
		_write(&labelLength, sizeof(uint32_t));
		_write(_header.c_str(), labelLength * sizeof(char));
	} else {
		uint32_t zero32 = 0;
		_write(&zero32, sizeof(uint32_t));
	}
};

void TGlfWriter::open(const std::string &Filename, const genometools::TChromosomes &Chrs, const std::string &Header){
	_filename = Filename;
	_gzfp     = nullptr;
	_gzfp     = gzopen(_filename.c_str(), "wb");

	if (_gzfp == nullptr) UERROR("Failed to open file '", _filename, "' for writing!");
	_index.clear();

	// write header
	_writeHeader(Header);
};

void TGlfWriter::newChromosome(const genometools::TChromosome &chromosome) {
	const uint8_t zero8 = 0;
	_write(&zero8, sizeof(uint8_t));

	// add chromosome info to index file
	_index.addChromosme(chromosome, _positionInFile);

	// write new chromosome: length of label, label, refId, length of ref sequence, ploidy
	const uint32_t labelLength = chromosome.name().size();
	_write(labelLength);
	_write(chromosome.name().c_str(), chromosome.name().length() * sizeof(char));
	_write(chromosome.refID());
	_write(chromosome.length());
	_write(chromosome.ploidy()); // TODO: I get an "uninitialized variable" error with valgrind. Why?

	// set oldPos
	_oldPos = 0;
};

void TGlfWriter::writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual,
			   const GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods) {
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
	if (_index.curChr().isHaploid()) {
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
	_write(glfValues.data(), _index.curChrNumLikelihoodValues() * sizeof(genometools::HighPrecisionPhredIntProbability));
}
void TGlfWriter::close() {
	if(_gzfp){
		_index.writeChromosmes(_filename);
	}
};

//---------------------------------
// TGlfReader
//---------------------------------

bool TGlfReader::_readChr() {
	// store position in file before read chromosome info
	auto posBeforeChr = _positionInFile - 1; // need to remove record type reading

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

	const auto refId  = _read<uint32_t>(); // currently not used
	const auto length = _read<uint32_t>();
	const auto ploidy = _read<uint8_t>();

	// check if chromosome is in index and set cur	
	if(_hasIndex){
		_index.setCurChromosomeToAndCheck(std::string_view(name, len), length, ploidy, posBeforeChr);
	} else {
		_index.addChromosme(std::string_view(name, len), length, ploidy, posBeforeChr);
	}

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
	      _index.curChrNumLikelihoodValues() * sizeof(genometools::HighPrecisionPhredIntProbability));
	_genotypeLikelihoodsGLF.type = _index.curChr().isHaploid() ? Ploidy::haploid : Ploidy::diploid;
};

void TGlfReader::setFilename(const std::string &Filename) { _filename = Filename; };

void TGlfReader::open(const std::string &Filename, bool HasIndex) {
	setFilename(Filename);
	_open(HasIndex);
};

void TGlfReader::_open(bool HasIndex) {	
	if (_gzfp) UERROR(_filename, " is already open!");
	_gzfp = gzopen(_filename.c_str(), "rb");

	if (_gzfp == nullptr) UERROR("Failed to open file '", _filename, "' for reading!");
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
	if(!_hasIndex){
		_index.readChromosomes(_filename + ".idx");
	}

	// read info of first chromosome
	_readRecordType();
	if (_recordType != 0)
		UERROR("GLF file does not start with chromosome entry. Are you using GLF files produced with an earlier version of ATLAS?");
	_readChr();
};

const genometools::TChromosomes& TGlfReader::chromosomes(){
	assert(_hasIndex);
	return _index.chromosomes();	
}

const genometools::TChromosome &TGlfReader::curChromosome() {
	assert(_hasIndex);
	return _index.curChr();
}

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

bool TGlfReader::jumpToChr(uint32_t RefID) {
	_index.seekChr(RefID);		
	if(gzseek64(_gzfp, _index.curChrPositionInFile(), SEEK_SET) < 0){ return false; }

	// read record type
	if (!_readRecordType()){ return false; }
	if(_recordType != 0){
		UERROR("File index (position in file) of chromosome '", _index.curChr().name(), "' does not point to a valid position within the GLF file!");
	}

	_readChr();
	return readNext();
};

bool TGlfReader::jumpToNextChr() {
	if(_index.curChrIsLast()) { return false; }
	_index.jumpToNextChromosome();
	jumpToChr(_index.curChr().refID());
};

bool TGlfReader::readNextWindow(std::vector<TGLFLikelihoods> &genoLikelihoods, uint32_t refId, uint32_t start,
				uint32_t end) {
	if (_eof) return false;

	// move to correct chromosome
	assert(_hasIndex);
	_index.seekChr(refId);
	gzseek64(_gzfp, _index.curChrPositionInFile(), SEEK_SET);
	if(!_readRecordType() || _recordType != 0){
		UERROR("File index (position in file) of chromosome with ref ID ", refId, " does not point to a valid position within the GLF file!");
	}
	_readChr();

	// jump to first position in window
	//  if there is a position 0 on the first chromosome parsed and windows starts right there
	//  -> we will not have read anything -> do this now
	while (_recordType == 0 || (_position < start && _index.curChr().refID() == refId)) {
		if (!readNext()) return false;
	}

	// have we passed window?
	if (_position >= end) return false; // no data

	// We are at first position in window with data

	// ensure size of container, fill with missing Data
	genoLikelihoods.resize(end - start, TGLFLikelihoods{});

	// Assumes that windows are read in order: no jumping back!
	if (refId < _index.curChr().refID()) { return false; }

	while (_position < end && _index.curChr().refID() == refId) {
		// fill in genotype likelihoods of current position
		genoLikelihoods[_position - start] = _genotypeLikelihoodsGLF;
		// read next record
		if (!readNext()) break;         // reached eof
	}
	return true;
};

// printing
void TGlfReader::printChr() {
	std::cout << "CHROMOSOME: '" << _index.curChr().name() << "' of length " << _index.curChr().length() << " and ploidy "
		  << (int)_index.curChr().ploidy() << "\n";
	
	// print header
	std::cout << "chr\tpos\tdepth\tRMS(MQ)";
	if(_index.curChr().isHaploid()){
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
	std::cout << _index.curChr().name() << "\t" << _position + 1 << "\t" << _depth << "\t" << _RMS_mappingQual;
	for (size_t i = 0; i < _index.curChrNumLikelihoodValues(); ++i) std::cout << "\t" << _genotypeLikelihoodsGLF.data()[i];
	std::cout << "\n";
};

void TGlfReader::printToEnd() { // For debugging
	std::cout << "GLF version is " << _version << "\n";	

	// first chromosome	
	printChr();

	// now parse file
	std::string oldChr = _index.curChr().name();
	while (readNext()) {
		if (oldChr != _index.curChr().name()) {
			printChr();
			oldChr = _index.curChr().name();
		}
		printSite();
	}
};

void TGLFPrinter::run() {
	TGlfReader reader(coretools::instances::parameters().get<std::string>("glf"));
	reader.printToEnd();
}

void TGLFIndexer::run(){
	// open GLF for reading, do not open index

}

}; // namespace GLF
