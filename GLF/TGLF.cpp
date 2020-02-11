/*
 * TGLF.cpp
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#include "TGLF.h"

//----------------------------------------------------
// TGlfConverter
//----------------------------------------------------
TGlfConverter::TGlfConverter(){
	maxVal = 65535;
	minLikelihood = pow(10.0, (double) maxVal / -1000.0);

	//fill map
	likelihoodMap.resize(maxVal+1);
	for(int i=0; i<=maxVal; i++){
		likelihoodMap[i] = pow(10.0, (double) i / -1000.0);
	}
};

uint16_t TGlfConverter::toGlfFormat(double scaledLikelihood){
	if(scaledLikelihood < minLikelihood){
		return maxVal;
	} else {
		return round(-1000.0 * log10(scaledLikelihood));
	}
};

uint16_t TGlfConverter::log10ToGlfFormat(double log10ScaledLikelihood){
	uint32_t tmp = round(-1000.0 * log10ScaledLikelihood);
	if(tmp > maxVal) return maxVal;
	else return tmp;
};

uint16_t TGlfConverter::phredToGlfFormat(uint8_t phred){
	return phred * 100;
};

double TGlfConverter::toScaledLikelihood(uint16_t glfValue){
	if(glfValue > maxVal) return minLikelihood;
	return likelihoodMap[glfValue];
};

uint8_t TGlfConverter::toPhred(uint16_t glfValue){
	uint16_t tmp = round(glfValue / 100.0);
	if(tmp > 255) return 255;
	return tmp;
};

double TGlfConverter::toLog10(uint16_t glfValue){
	return glfValue / -1000.0;
};

double TGlfConverter::toLog(uint16_t glfValue){
	return glfValue * logOf10DividedByMinus1000;
};

//---------------------------------
//TGlfWriter
//---------------------------------
//PRIVATE
void TGlfWriter::init(){
	glfValues = new uint16_t[curChr.maxNumLikelihoodValues];
	oldPos = 0;
	recordType1 = one8 << 4;
};

void TGlfWriter::writeHeader(){
	write(version.c_str(), 4*sizeof(char));

	if(header.length() > 0){
		uint32_t labelLength = header.size();
		write(&labelLength, sizeof(uint32_t));
		write(header.c_str(), labelLength * sizeof(char));
	} else
		write(&zero32, sizeof(uint32_t));
};

//PUBLIC
void TGlfWriter::open(std::string Filename, std::string Header){
	filename = Filename;
	gzfp = NULL;
	gzfp = gzopen(filename.c_str(),"wb");

	if(gzfp == NULL)
		throw "Failed to open file '" + filename + "' for writing!";
	isOpen = true;
	curChr.clear();

	//write header
	header = Header;
	writeHeader();
};

void TGlfWriter::newChromosome(const std::string name, const uint16_t refId, const uint32_t length, const uint8_t ploidy){
	if(curChr.name != "")
		write(&zero8, sizeof(uint8_t));

	//save cur info
	curChr.update(name, refId, length, ploidy);

	//write new chromosome: length of label, label, refId, length of ref sequence, ploidy
	uint32_t labelLength = curChr.name.size();

	write(&labelLength, sizeof(uint32_t));
	write(curChr.name.c_str(), name.length() * sizeof(char));
	write(&curChr.refId, sizeof(uint32_t));
	write(&curChr.length, sizeof(uint32_t));
	write(&curChr.ploidy, sizeof(uint8_t)); // I get an "uninitialized varable" error with valgrind. Why?

	//set oldPos and curChr
	oldPos = 0;
};

void TGlfWriter::writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual, double* genotypeLikelihoods){
	//record type
	//TODO: add reference?
	write(&recordType1, sizeof(uint8_t));

	//offset
	offset = pos - oldPos;
	oldPos = pos;
	write(&offset, sizeof(uint32_t));

	//calculate likelihoods in GLF format
	//Note: genotype likelihoods are given for the 10 diploid genotypes!!
	//TODO: maybe do in GLFChromosomes?
	if(curChr.ploidy == 1){
		double maxLL = genotypeLikelihoods[AA];
		if(genotypeLikelihoods[CC] > maxLL) maxLL = genotypeLikelihoods[CC];
		if(genotypeLikelihoods[GG] > maxLL) maxLL = genotypeLikelihoods[GG];
		if(genotypeLikelihoods[TT] > maxLL) maxLL = genotypeLikelihoods[TT];

		//normalize and scale to uint16
		glfValues[0] = converter.toGlfFormat(genotypeLikelihoods[AA] / maxLL);
		glfValues[1] = converter.toGlfFormat(genotypeLikelihoods[CC] / maxLL);
		glfValues[2] = converter.toGlfFormat(genotypeLikelihoods[GG] / maxLL);
		glfValues[3] = converter.toGlfFormat(genotypeLikelihoods[TT] / maxLL);
	} else {
		//maxLL
		double maxLL = genotypeLikelihoods[0];
		for(int i=1; i<curChr.numLikelihoodValues; ++i){
			if(genotypeLikelihoods[i] > maxLL)
				maxLL = genotypeLikelihoods[i];
		}

		//normalize and scale to uint16
		for(int i=0; i<curChr.numLikelihoodValues; ++i){
			glfValues[i] = converter.toGlfFormat(genotypeLikelihoods[i] / maxLL);
		}
	}

	//write maxLL as uint16_t
	//uint16_t maxLL_int = converter.toGlfFormat(maxLL);
	//write(&maxLL_int, sizeof(uint16_t));

	//write depth as uint16_t
	if(depth > 65535) depth = 65535;
	uint16_t tmp = depth;
	write(&tmp, sizeof(uint16_t));

	//root mean square of mapping qualities
	write(&RMS_mappingQual, sizeof(uint8_t));

	//genotype likelihoods
	write(glfValues, curChr.numLikelihoodValues*sizeof(uint16_t));
};

//---------------------------------
//TGlfReader
//---------------------------------
//PRIVATE
void TGlfReader::init(){
	isOpen = false;
	reachedEndOfChr = true;
	recordType = 99;
	offset = 0;
	_position = 0;
	//maxLL = 0;
	_depth = 0;
	_RMS_mappingQual = 0;
	positionInFile = 0;
	tmpInt8 = 0;
	_lenRead = 0;
	_eof = true;

	_genotypeLikelihoodsGLF_missingData = new uint16_t[curChr.maxNumLikelihoodValues];
	for(int i=0; i<curChr.maxNumLikelihoodValues; ++i)
		_genotypeLikelihoodsGLF_missingData[i] = 0;
};

bool TGlfReader::fillPointerToChr(uint32_t refId, TGlfChromosome* & chr){
	if(curChr.refId == refId){
		chr = &curChr;
		return true;
	} else {
		std::map< uint32_t, TGlfChromosome >::iterator it = _chromosomesAlreadyParsed.find(refId);
		if(it == _chromosomesAlreadyParsed.end()){
			chr = nullptr;
			return false;
		} else {
			chr = &it->second;
			return true;
		}
	}
};

bool TGlfReader::readChr(){
	//store current chromosome name in list of chromosomes parsed
	if(curChr.name != "")
		_chromosomesAlreadyParsed.emplace(curChr.refId, curChr);

	//read chromosome info
	//name
	uint32_t len;
	if(!read(&len, sizeof(uint32_t))){
		_eof = true;
		return false;
	}
	char* tmp = new char[len];
	read(tmp, len*sizeof(char));
	std::string name(tmp, len);
	delete[] tmp;

	//refId
	uint16_t refId;
	read(&refId, sizeof(uint32_t));

	//length
	uint32_t length;
	read(&length, sizeof(uint32_t));

	//ploidy
	uint8_t ploidy;
	read(&ploidy, sizeof(uint8_t));

	//update
	curChr.update(name, refId, length, ploidy);

	//set first position = 0
	_position = 0;

	return true;
};

bool TGlfReader::readRecordType(){
	if(!read(&tmpInt8, sizeof(uint8_t))){
		_eof = true;
		return false;
	}
	recordType = tmpInt8 >> 4;
	if(recordType > 1) throw "Unknown record type in file '" + filename + "'!";
	return true;
};

void TGlfReader::readSNPRecord(){
	//read data of a single position
	//offset
	read(&offset, sizeof(uint32_t));
	_position += offset;

	//maxLL and depth (both uint16)
	//read(&maxLL, sizeof(uint16_t));
	read(&_depth, sizeof(uint16_t));

	//root mean square of mapping qualities
	read(&_RMS_mappingQual, sizeof(uint8_t));
	//RMS_mappingQual = (int) tmpInt8; DO WE NEED THIS??

	//genotype likelihoods
	read(_genotypeLikelihoodsGLF, curChr.numLikelihoodValues*sizeof(uint16_t));
};


//PUBLIC
void TGlfReader::setFilename(std::string Filename){
	filename = Filename;
};

void TGlfReader::open(std::string Filename){
	setFilename(Filename);
	open();
};

void TGlfReader::open(){
	gzfp = NULL;
	gzfp = gzopen(filename.c_str(),"rb");

	if(gzfp == NULL)
		throw "Failed to open file '" + filename + "' for reading!";
	isOpen = true;
	curChr.clear();
	positionInFile = 0;

	//parse header
	//version
	char buffer[4];
	read(buffer, 4*sizeof(char));
	std::string fileVersion;
	fileVersion.assign(buffer, 4);

	if(fileVersion != version)
		throw "Non-supported GLF version '" + fileVersion + "! Currently only version '" + version + "' is supported.";
	version = fileVersion;

	read(&HeaderLen, sizeof(uint32_t));

	if(HeaderLen > 0){
		char* header = new char[HeaderLen];
		read(&header, HeaderLen*sizeof(char));
	}
	_eof = false;

	//read info of first chromosome
	_chromosomesAlreadyParsed.clear();
	readChr();
};

void TGlfReader::rewind(){
	//go back to beginning of file
	close();
	open();
};

bool TGlfReader::readNext(){
	//read record type
	if(!readRecordType()) return false;
	if(recordType == 0){
		readChr();
		return readNext();
	} else if(recordType == 1){
		readSNPRecord();
		return true;
	} else throw "Unknown record type in file '" + filename + "'!";
};

bool TGlfReader::jumpToEndOfChr(){
	//read record type
	if(!readRecordType()) return false;

	//tmp variables
	while(recordType != 0){
		//skipRecord();
		readSNPRecord();

		if(!readRecordType()) return false;
	}

	return true;
};

bool TGlfReader::jumpToNextChr(){
	if(!jumpToEndOfChr()){
		return false;
	}
	readChr();
	return readNext();
};

bool TGlfReader::readNextWindow(std::vector<uint16_t*> & genoLikelihoods, const uint32_t refId, const uint32_t start, const uint32_t end){
	//Assumes that windows are read in order: no jumping back!
	if(_eof) return false;
	if(refId < curChr.refId) return false;

	//move to correct chromosome
	while(curChr.refId < refId){
		jumpToNextChr();
	}

	//jump to first position in window
	while(_position < start && curChr.refId == refId){
		readNext();
	}

	//have we passed window?
	if(_position >= end)
		return false; //no data

	//We are at first position in window with data
	uint32_t i = start;
	int index = 0;

	while(_position < end && curChr.refId == refId){
		//fill in missing positions before
		for(; i<_position; ++i, ++index)
			memcpy(genoLikelihoods[index], _genotypeLikelihoodsGLF_missingData, curChr.numLikelihoodValues*sizeof(uint16_t));

		memcpy(genoLikelihoods[index], _genotypeLikelihoodsGLF, curChr.numLikelihoodValues*sizeof(uint16_t));
		++index; ++i;

		//read next record
		readNext();
		if(curChr.refId != refId){
			for(; i<end; ++i, ++index)
				memcpy(genoLikelihoods[index], _genotypeLikelihoodsGLF_missingData, curChr.numLikelihoodValues*sizeof(uint16_t));
		}
	}

	return true;
};

void TGlfReader::fillGenotypeLikelihoodsGLF(uint16_t* destination){
	//assumes pointer points to
	memcpy(destination, _genotypeLikelihoodsGLF, curChr.numLikelihoodValues*sizeof(uint16_t));
};

void TGlfReader::fillGenotypeLikelihoods(double* destination, TGlfConverter* converter){
	for(int i=0; i<curChr.numLikelihoodValues; ++i){
		destination[i] = converter->toScaledLikelihood(_genotypeLikelihoodsGLF[i]);
	}
};

//printing
void TGlfReader::printChr(){
	std::cout << "CHROMOSOME: '" << curChr.name << "' of length " << curChr.length << " and ploidy " << (int) curChr.ploidy << "\n";
};

void TGlfReader::printSite(){
	//std::cout << curChr.name << "\t" << position << "\t" << maxLL << "\t" << depth << "\t" << RMS_mappingQual;
	std::cout << curChr.name << "\t" << _position << "\t" << _depth << "\t" << _RMS_mappingQual;
	for(int i=0; i<curChr.numLikelihoodValues; ++i)
		std::cout << "\t" << unsigned(_genotypeLikelihoodsGLF[i]);
	std::cout << "\n";
};

void TGlfReader::printToEnd(){ //For debugging
	//first print header
	std::cout << version << "\n";
	std::cout << header << "\n";
	printChr();

	//now parse file
	std::string oldChr = curChr.name;
	while(readNext()){
		if(oldChr != curChr.name){
			printChr();
			oldChr = curChr.name;
		}
		printSite();
	}
};
