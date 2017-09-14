/*
 * TGLF.cpp
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#include "TGLF.h"

//---------------------------------
//TGlfWriter
//---------------------------------
//PRIVATE
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
	positionInFile = 0;
	curChr = "";

	//write header
	header = Header;
	writeHeader();
};

void TGlfWriter::newChromosome(std::string name, uint32_t length){
	if(curChr != "")
		write(&zero8, sizeof(uint8_t));

	//write new chromosome: length of label, label, length of ref sequence
	uint32_t labelLength = name.size();

	write(&labelLength, sizeof(uint32_t));
	write(name.c_str(), name.length() * sizeof(char));
	write(&length, sizeof(uint32_t));

	//set oldPos and curChr
	oldPos = 0;
	curChr = name;
};

void TGlfWriter::writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual, uint8_t* genoQualities, uint32_t & maxLL){
	//record type
	//TODO: add reference
	write(&recordType1, sizeof(uint8_t));

	//offset
	offset = pos - oldPos;
	oldPos = pos;
	write(&offset, sizeof(uint32_t));

	//maxLL (capped at 255) and depth as one uint32_t
	if(maxLL > 255) maxLL = 255;
	uint32_t tmp = maxLL << 24;
	tmp = tmp | depth;
	write(&tmp, sizeof(uint32_t));

	//root mean square of mapping qualities
	write(&RMS_mappingQual, sizeof(uint8_t));

	//genotype likelihoods
	//TODO: test if using more accuracy on qualities matters
	write(genoQualities, 10*sizeof(uint8_t));
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
	position = 0;
	maxLL = 0;
	depth = 0;
	RMS_mappingQual = 0;
	curChr = "";
	_chrLength = 0;
	depth_mask = 0xFFFFFF;
	tmpInt32 = 0;
	tmpInt8 = 0;
	_lenRead = 0;
	_i = 0;
	_eof = true;
	genotypeQualitiesMissingData = new int[10];
	for(int i=0; i<10; ++i)
		genotypeQualitiesMissingData[i] = 0;
};

bool TGlfReader::readChr(){
	//store current chromosome name in list of chromosomes parsed
	if(curChr != "")
		chromosomesAlreadyParsed.push_back(curChr);

	//read chromosome info
	uint32_t len;
	if(!read(&len, sizeof(uint32_t))){
		_eof = true;
		return false;
	}
	char* tmp = new char[len];
	read(tmp, len*sizeof(char));
	curChr.assign(tmp, len);

	read(&tmpInt32, sizeof(uint32_t));
	_chrLength = tmpInt32;
	position = 0;

	return true;
};

bool TGlfReader::chromosomeParsed(std::string & chr){
	for(std::vector<std::string>::iterator it=chromosomesAlreadyParsed.begin(); it!=chromosomesAlreadyParsed.end(); ++it){
		if(*it == chr)
			return true;
	}
	return false;
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
	position += offset;

	//maxLL and depth
	read(&tmpInt32, sizeof(uint32_t));
	depth = tmpInt32 & depth_mask;
	maxLL = tmpInt32 >> 24;

	//root mean square of mapping qualities
	read(&tmpInt8, sizeof(uint8_t));
	RMS_mappingQual = (int) tmpInt8;

	//genotype likelihoods
	read(tmpInt8_10, 10*sizeof(uint8_t));

	for(_i=0; _i<10; ++_i){
		genotypeQualities[_i] = (int) tmpInt8_10[_i];
	}
};


//PUBLIC
void TGlfReader::open(std::string Filename){
	filename = Filename;
	gzfp = NULL;
	gzfp = gzopen(filename.c_str(),"rb");

	if(gzfp == NULL)
		throw "Failed to open file '" + filename + "' for reading!";
	isOpen = true;
	positionInFile = 0;
	curChr = "";

	//parse header
	//version
	char buffer[4];
	read(buffer, 4*sizeof(char));
	version = buffer;

	read(&HeaderLen, sizeof(uint32_t));

	if(HeaderLen > 0){
		char* header = new char[HeaderLen];
		read(&header, HeaderLen*sizeof(char));
	}
	_eof = false;

	//read info of first chromosome
	chromosomesAlreadyParsed.clear();
	readChr();
};

void TGlfReader::rewind(){
	//go back to beginning of file
	close();
	open(filename);
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
		skipRecord();
		if(!readRecordType()) return false;
	}
	return true;
};

bool TGlfReader::readNextWindow(int** genoLikelihoods, std::string chr, long start, long end){
	//Assumes that windows are read in order: no jumping back!
	if(_eof) return false;
	if(chromosomeParsed(chr)) return false;

	//move to correct chromosome
	if(curChr != chr){
		while(curChr != chr){
			jumpToEndOfChr();
			readChr();
		}
		if(!readRecordType()) return false;
		if(recordType == 0){
			//means no data on this chromosome
			readChr();
			return false;
		} else readSNPRecord();
	}

	//jump to first position in window
	while(position < start){
		if(!readRecordType()) return false;
		if(recordType == 0){
			//means no data on this chromosome
			readChr();
			return false;
		} else readSNPRecord();
	}

	//have we passed window?
	if(position >= end)
		return false; //no data

	//We are at first position in window with data
	long i = start;
	int index = 0;

	while(position < end){
		//fill in missing positions before
		for(; i<position; ++i, ++index)
			memcpy(genoLikelihoods[index], genotypeQualitiesMissingData, 10*sizeof(int));

		//now add data
		//std::cout << "i = " << i << "\tindex = " << index << "\tposition = " << position << "\tgenotypeQualities[0] = " << genotypeQualities[0] << std::endl;

		memcpy(genoLikelihoods[index], genotypeQualities, 10*sizeof(int));
		++index; ++i;

		//read next record
		if(!readRecordType() || recordType == 0){
			for(; i<end; ++i, ++index)
				memcpy(genoLikelihoods[index], genotypeQualitiesMissingData, 10*sizeof(int));
			readChr();
			break;
		}
		readSNPRecord();
	}

	return true;
}

//printing
void TGlfReader::printChr(){
	std::cout << "CHROMOSOME: '" << curChr << "' of length " << _chrLength << "\n";
};

void TGlfReader::printSite(){
	std::cout << curChr << "\t" << position << "\t" << maxLL << "\t" << depth << "\t" << RMS_mappingQual;
	for(int i=0; i<10; ++i)
		std::cout << "\t" << genotypeQualities[i];
	std::cout << "\n";
};

void TGlfReader::printToEnd(){ //For debugging
	//first print header
	std::cout << version << "\n";
	std::cout << header << "\n";
	printChr();

	//now parse file
	std::string oldChr = curChr;
	while(readNext()){
		if(oldChr != curChr){
			printChr();
			oldChr = curChr;
		}
		printSite();
	}
};
