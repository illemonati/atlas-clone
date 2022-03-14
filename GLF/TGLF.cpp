/*
 * TGLF.cpp
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#include "TGLF.h"
#include "GenotypeTypes.h"
#include "debugtools.h"

namespace GLF{
using namespace GenotypeLikelihoods;

//----------------------------------------------------
//TGlfChromosome
//----------------------------------------------------
TGlfChromosome::TGlfChromosome(){
	_refId = 0;
	_length = 0;
	_isHaploid = false;
	_numLikelihoodValues = 10;
	_name = "";
};

TGlfChromosome::TGlfChromosome(const std::string & Name, uint32_t Length, const uint8_t & Ploidy){
	_name = Name;
	_length = Length;
	_refId = 0;
	_setPloidy(Ploidy);
};

void TGlfChromosome::_setPloidy(const uint8_t & Ploidy){
	if(Ploidy < 1 || Ploidy > 2)
		throw "Currently GLFs only support ploidies 1 and 2 (not " + coretools::str::toString((int) Ploidy) + ")!";
	if(Ploidy == 1){
		_isHaploid = true;
		_numLikelihoodValues = 4;
	} else {
		_isHaploid = false;
		_numLikelihoodValues = 10;
	}
};

void TGlfChromosome::update(const std::string & Name, uint16_t RefId, uint32_t Length, const uint8_t & Ploidy){
	_name = Name;
	_refId = RefId;
	_length = Length;
	_setPloidy(Ploidy);
};


void TGlfChromosome::clear(){
	_name = "";
	_refId = 0;
	_length = 0;
	_isHaploid = false;
};


//---------------------------------
//TGlfWriter
//---------------------------------
void TGlfWriter::_init(){
    _oldPos = 0;
    _recordType1 = _one8 << 4;
};

void TGlfWriter::_writeHeader(){
    _write(_version.c_str(), 4 * sizeof(char));

	if(_header.length() > 0){
		uint32_t labelLength = _header.size();
        _write(&labelLength, sizeof(uint32_t));
        _write(_header.c_str(), labelLength * sizeof(char));
	} else
        _write(&_zero32, sizeof(uint32_t));
};

void TGlfWriter::open(const std::string& Filename){
	open(Filename, "");
};

void TGlfWriter::open(const std::string& Filename, const std::string& Header){
    _filename = Filename;
    _gzfp = nullptr;
    _gzfp = gzopen(_filename.c_str(), "wb");

	if(_gzfp == nullptr)
		throw "Failed to open file '" + _filename + "' for writing!";
    _isOpen = true;
	_curChr.clear();

	//write header
	_header = Header;
    _writeHeader();
};

void TGlfWriter::newChromosome(const BAM::TChromosome & chromosome){
    _write(&_zero8, sizeof(uint8_t));

	//save cur info
	_curChr.update(chromosome.name, chromosome.refID(), chromosome.length, chromosome.ploidy);

	//write new chromosome: length of label, label, refId, length of ref sequence, ploidy
	uint32_t labelLength = _curChr.name().size();
    _write(labelLength);
    _write(_curChr.name().c_str(), _curChr.name().length() * sizeof(char));
    _write(_curChr.refId());
    _write(_curChr.length());
    _write(_curChr.ploidy()); // TODO: I get an "uninitialized variable" error with valgrind. Why?

	//set oldPos and curChr
	_oldPos = 0;
};

void TGlfWriter::writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods){
	using genometools::Genotype;
	//record type
	//TODO: add reference?
    _write(&_recordType1, sizeof(uint8_t));

	//offset
	_offset = pos - _oldPos;
    _oldPos = pos;
    _write(&_offset, sizeof(uint32_t));

	//calculate likelihoods in GLF format
	//Note: genotype likelihoods are given for the 10 diploid genotypes!!
	if(_curChr.isHaploid()){
		coretools::Probability maxLik = genotypeLikelihoods[Genotype::AA];
		if(genotypeLikelihoods[Genotype::CC] > maxLik) maxLik = genotypeLikelihoods[Genotype::CC];
		if(genotypeLikelihoods[Genotype::GG] > maxLik) maxLik = genotypeLikelihoods[Genotype::GG];
		if(genotypeLikelihoods[Genotype::TT] > maxLik) maxLik = genotypeLikelihoods[Genotype::TT];

		//normalize and scale to uint16
		_glfValues[0] = genotypeLikelihoods[Genotype::AA] / maxLik;
        _glfValues[1] = genotypeLikelihoods[Genotype::CC] / maxLik;
        _glfValues[2] = genotypeLikelihoods[Genotype::GG] / maxLik;
        _glfValues[3] = genotypeLikelihoods[Genotype::TT] / maxLik;
	} else {
		//ploidy is 2
		coretools::Probability maxLik = genotypeLikelihoods.max();

		//normalize and scale to genometools::HighPrecisionPhredIntProbability
		for(auto g = Genotype::min; g < Genotype::max; ++g){
			_glfValues[genometools::index(g)] = genotypeLikelihoods[g] / maxLik;
		}
	}

	//write maxLL as uint16_t
	//uint16_t maxLL_int = converter.toGlfFormat(maxLL);
	//write(&maxLL_int, sizeof(uint16_t));

	//write depth as uint16_t
	if(depth > 65535) depth = 65535;
	uint16_t tmp = depth;
    _write(&tmp, sizeof(uint16_t));

	//root mean square of mapping qualities
    _write(&RMS_mappingQual, sizeof(uint8_t));

	//genotype likelihoods
    _write(_glfValues.data(), _curChr.numLikelihoodValues() * sizeof(genometools::HighPrecisionPhredIntProbability));
};

//---------------------------------
//TGlfReader
//---------------------------------
void TGlfReader::_init(){
    _isOpen = false;
    _reachedEndOfChr = true;
    _recordType = 99;
	_offset = 0;
	_position = 0;
	//maxLL = 0;
	_depth = 0;
	_RMS_mappingQual = 0;
    _positionInFile = 0;
    _tmpInt8 = 0;
	_lenRead = 0;
	_eof = true;

	for(int i=0; i < 10; ++i){
		_genotypeLikelihoodsGLF_missingData[i] = genometools::HighPrecisionPhredIntProbability::highest();
	}
};

TGlfChromosome* TGlfReader::pointerToChr(uint32_t refId){
	if(_curChr.refId() == refId){
		return &_curChr;
	} else {
		auto it = _chromosomesAlreadyParsed.find(refId);
		if(it == _chromosomesAlreadyParsed.end()){
			return nullptr;
		} else {
			return &it->second;
		}
	}
};

bool TGlfReader::fillPointerToChr(uint32_t refId, TGlfChromosome* & chr){
	chr = pointerToChr(refId);
	if(chr == nullptr) return false;
	else return true;
};

bool TGlfReader::_readChr(){
	//store current chromosome name in list of chromosomes parsed
	if(!_curChr.name().empty())
		_chromosomesAlreadyParsed.emplace(_curChr.refId(), _curChr);

	//read chromosome info
	//name
	uint32_t len;
	if(!_read(&len, sizeof(uint32_t))){
		_eof = true;
		return false;
	}
	char* tmp = new char[len];
    _read(tmp, len * sizeof(char));
	std::string name(tmp, len);
	delete[] tmp;

	//refId
	uint16_t refId;
    _read(&refId, sizeof(uint32_t));

	//length
	uint32_t length;
    _read(&length, sizeof(uint32_t));

	//ploidy
	uint8_t ploidy;
    _read(&ploidy, sizeof(uint8_t));

	//update
	_curChr.update(name, refId, length, ploidy);

	//set first position = 0
	_position = 0;

	return true;
};

bool TGlfReader::_readRecordType(){
	if(!_read(&_tmpInt8, sizeof(uint8_t))){
		_eof = true;
		return false;
	}
    _recordType = _tmpInt8 >> 4;
	if(_recordType > 1) throw "Unknown record type in file '" + _filename + "'!";
	return true;
};

void TGlfReader::_readSNPRecord(){
	//read data of a single position
	//offset
    _read(&_offset, sizeof(uint32_t));
	_position += _offset;

	//maxLL and depth (both uint16)
	//read(&maxLL, sizeof(uint16_t));
    _read(&_depth, sizeof(uint16_t));

	//root mean square of mapping qualities
    _read(&_RMS_mappingQual, sizeof(uint8_t));
	//RMS_mappingQual = (int) tmpInt8; DO WE NEED THIS??

	//genotype likelihoods
    _read(_genotypeLikelihoodsGLF.data(), _curChr.numLikelihoodValues() * sizeof(genometools::HighPrecisionPhredIntProbability));
};

void TGlfReader::setFilename(const std::string& Filename){
    _filename = Filename;
};

void TGlfReader::open(const std::string& Filename){
	setFilename(Filename);
    _open();
};

void TGlfReader::_open(){
    _gzfp = nullptr;
    _gzfp = gzopen(_filename.c_str(), "rb");

	if(_gzfp == nullptr)
		throw "Failed to open file '" + _filename + "' for reading!";
    _isOpen = true;
	_curChr.clear();
    _positionInFile = 0;

    //parse header
	//version
	char buffer[4];
    _read(buffer, 4 * sizeof(char));
	std::string fileVersion;
	fileVersion.assign(buffer, 4);

	if(fileVersion != _version)
		throw "Non-supported GLF version '" + fileVersion + "! Currently only version '" + _version + "' is supported.";
    _version = fileVersion;

    _read(&_HeaderLen, sizeof(uint32_t));

	if(_HeaderLen > 0){
		char* header = new char[_HeaderLen];
        _read(&header, _HeaderLen * sizeof(char));
	}
	_eof = false;

	//read info of first chromosome
	_chromosomesAlreadyParsed.clear();
	_readRecordType();
	if (_recordType != 0)
	    throw "GLF file does not start with chromosome entry. The GLF format has changed with ATLAS 1.0, are you using GLF files produced with an earlier version?";
    _readChr();
};

void TGlfReader::rewind(){
	//go back to beginning of file
	close();
    _open();
};

bool TGlfReader::readNext(){
	//read record type
	if(!_readRecordType()) return false;
	if(_recordType == 0){
        _readChr();
		return readNext();
	} else if(_recordType == 1){
        _readSNPRecord();
		return true;
	} else throw "Unknown record type in file '" + _filename + "'!";
};

bool TGlfReader::_jumpToEndOfChr(){
	//read record type
	if(!_readRecordType()) return false;

	//tmp variables
	while(_recordType != 0){
		//skipRecord();
        _readSNPRecord();

		if(!_readRecordType()) return false;
	}

	return true;
};

bool TGlfReader::jumpToNextChr(){
	if(!_jumpToEndOfChr()){
		return false;
	}
    _readChr();
	return readNext();
};

bool TGlfReader::readNextWindow(std::vector<GLFLikelihoods> & genoLikelihoods, uint32_t refId, uint32_t start, uint32_t end){
	if(_eof) return false;

	//move to correct chromosome
	while(_curChr.refId() < refId){
        if (!jumpToNextChr())
            return false; // chromosome id doesn't exist (reached eof)
    }

	//jump to first position in window
    // if there is a position 0 on the first chromosome parsed and windows starts right there
    // -> we will not have read anything -> do this now
	while(_recordType == 0 || (_position < start && _curChr.refId() == refId)){
		if (!readNext())
            return false;
    }

	//have we passed window?
	if(_position >= end)
        return false; //no data

	//We are at first position in window with data
	uint32_t i = start;
	int index = 0;

	//ensure size of container
	genoLikelihoods.resize(end - start + 1);

    //Assumes that windows are read in order: no jumping back!
    if(refId < _curChr.refId()){
        return false;
    }

	while(_position < end && _curChr.refId() == refId) {
        //fill in missing positions before
        for (; i < _position; ++i, ++index){
        	genoLikelihoods[index] = _genotypeLikelihoodsGLF_missingData;
        }

        // fill in genotype likelihoods of current position
        genoLikelihoods[index] = _genotypeLikelihoodsGLF;
		++index; ++i;

		//read next record
		if (!readNext())
		    break; // reached eof
        if(_curChr.refId() != refId){ // reached next chromosome
        	for(; index < genoLikelihoods.size(); ++index){
            	genoLikelihoods[index] = _genotypeLikelihoodsGLF_missingData;
            }
		}
	}

	// fill sites that were in the end of the window
    for(; index < genoLikelihoods.size(); ++index){
    	genoLikelihoods[index] = _genotypeLikelihoodsGLF_missingData;
    }

	return true;
};

//printing
void TGlfReader::printChr(){
	std::cout << "CHROMOSOME: '" << _curChr.name() << "' of length " << _curChr.length() << " and ploidy " << (int) _curChr.ploidy() << "\n";
};

void TGlfReader::printSite(){
	//std::cout << curChr.name << "\t" << position << "\t" << maxLL << "\t" << depth << "\t" << RMS_mappingQual;
	// print position as 1-based, internally it is 0-based
	std::cout << _curChr.name() << "\t" << _position + 1 << "\t" << _depth << "\t" << _RMS_mappingQual;
	for(int i=0; i < _curChr.numLikelihoodValues(); ++i)
		std::cout << "\t" << _genotypeLikelihoodsGLF[i];
	std::cout << "\n";
};

void TGlfReader::printToEnd(){ //For debugging
	//first print header
	std::cout << _version << "\n";
	std::cout << _header << "\n";
	printChr();

	//now parse file
	std::string oldChr = _curChr.name();
	while(readNext()){
		if(oldChr != _curChr.name()){
			printChr();
			oldChr = _curChr.name();
		}
		printSite();
	}
};

}; //end namespace
