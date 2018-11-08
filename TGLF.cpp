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
	++curChrNumber;
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
	curChrNumber = 0;
	positionInFile = 0;
	_chrLength = 0;
	depth_mask = 0xFFFFFF;
	tmpInt32 = 0;
	tmpInt8 = 0;
	_lenRead = 0;
	_eof = true;

	for(int i=0; i<10; ++i)
		genotypeQualitiesMissingData[i] = 0;
};

std::string TGlfReader::getNameOfParsedChr(int chrNumber){
	if(curChrNumber == chrNumber)
		return curChr;
	else if(curChrNumber > chrNumber)
		return chromosomesAlreadyParsed[chrNumber].first;
	else throw "TGlfReader does not know name of chromosome " + toString(chrNumber) + ": chromosome not yet read!";
};

long TGlfReader::getLengthOfParsedChr(int chrNumber){
	if(curChrNumber == chrNumber)
		return _chrLength;
	else if(curChrNumber > chrNumber)
		return chromosomesAlreadyParsed[chrNumber].second;
	else throw "TGlfReader does not know length of chromosome " + toString(chrNumber) + ": chromosome not yet read!";
};

bool TGlfReader::readChr(){
	//store current chromosome name in list of chromosomes parsed
	if(curChr != "")
		chromosomesAlreadyParsed.push_back(std::pair<std::string, long>(curChr, _chrLength));

	//read chromosome info
	uint32_t len;
	if(!read(&len, sizeof(uint32_t))){
		_eof = true;
		return false;
	}
	char* tmp = new char[len];
	read(tmp, len*sizeof(char));
	curChr.assign(tmp, len);
	++curChrNumber;

	read(&tmpInt32, sizeof(uint32_t));
	_chrLength = tmpInt32;
	position = 0;

	delete[] tmp;
	return true;
};

bool TGlfReader::chromosomeParsed(std::string & chr){
	for(std::vector< std::pair<std::string, long> >::iterator it=chromosomesAlreadyParsed.begin(); it!=chromosomesAlreadyParsed.end(); ++it){
		if(it->first == chr)
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
	read(genotypeQualities, 10*sizeof(uint8_t));
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
	curChr = "";
	curChrNumber = 0;
	positionInFile = 0;

	//parse header
	//version
	char buffer[4];
	read(buffer, 4*sizeof(char));
	version.assign(buffer, 4);

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

bool TGlfReader::readNextWindow(std::vector<uint8_t*> & genoLikelihoods, std::string chr, long start, long end){
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
};

void TGlfReader::fillGenotypeQualities(uint8_t* destination){
	//assumes pointer points to
	memcpy(destination, genotypeQualities, 10*sizeof(uint8_t));
};

//printing
void TGlfReader::printChr(){
	std::cout << "CHROMOSOME: '" << curChr << "' of length " << _chrLength << "\n";
};

void TGlfReader::printSite(){
	std::cout << curChr << "\t" << position << "\t" << maxLL << "\t" << depth << "\t" << RMS_mappingQual;
	for(int i=0; i<10; ++i)
		std::cout << "\t" << unsigned(genotypeQualities[i]);
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

//----------------------------------------------------
//TGlfMultiReader
//----------------------------------------------------
TGlfMultiReader::TGlfMultiReader(){
	init();
};

TGlfMultiReader::TGlfMultiReader(std::vector<std::string> FileNames, TLog* logfile){
	init();
	openGLFs(FileNames, logfile);
};

TGlfMultiReader::TGlfMultiReader(TParameters & params, TLog* logfile){
	init();
	openGLFs(params, logfile);
};

void TGlfMultiReader::init(){
	numGLFs = 0;
	readersOpened = false;
	GLFs = NULL;
	data = NULL;
	GLFIsActive = NULL;
	hasData = NULL;
	dataInitialized = false;
	numActiveFiles = 0;

	_position = -1;
	_curChrNumber = 0;
	_curChrLength = 0;
	_curChrName = "";
	_numActiveFilesWithData = 0;

	for(int i=0; i<10; ++i)
			genotypeQualitiesMissingData[i] = 0;

};

TGlfMultiReader::~TGlfMultiReader(){
	if(dataInitialized){
		delete[] data;
		delete[] hasData;
	}

	if(readersOpened){
		delete[] GLFIsActive;
		closeGLF();
	}
};

void TGlfMultiReader::_openGLFs(TLog* logfile){
	numGLFs = GLFNames.size();
	GLFIsActive = new bool[numGLFs];

	//open files
	GLFs = new TGlfReader[numGLFs];
	readersOpened = true;
	logfile->startIndent("Opening GLF files:");
	int g = 0;
	for(std::vector<std::string>::iterator it=GLFNames.begin(); it != GLFNames.end(); ++it, ++g){
		logfile->listFlush("Opening GLF '" + *it + "' ...");
		GLFs[g].open(*it);
		logfile->done();
	}
	logfile->endIndent();

	_setAllInactive();
};


void TGlfMultiReader::openGLFs(const std::vector<std::string> & FileNames, TLog* logfile){
	GLFNames = FileNames;
	_openGLFs(logfile);
};

void TGlfMultiReader::openGLFs(TParameters & params, TLog* logfile){
	params.fillParameterIntoVector("glf", GLFNames, ',');
	_openGLFs(logfile);
};

void TGlfMultiReader::closeGLF(){
	if(readersOpened){
		//close all glf handlers
		for(int g=0; g<numGLFs; ++g)
			GLFs[g].close();

		delete[] GLFs;
		GLFNames.clear();
		numGLFs = 0;
		readersOpened = false;
	}
};

//-------------------------------------
//set active / inactive
//-------------------------------------
int TGlfMultiReader::_getGLFIndexFromName(const std::string & name){
	int index = 0;
	for(std::vector<std::string>::iterator it=GLFNames.begin(); it!=GLFNames.end(); ++it, ++index){
		if(*it == name) return index;
	}
	throw "GLF with name '" + name + "' not in TGlfMultiReader!";
};

void TGlfMultiReader::_setActive(const int index){
	if(index >= numGLFs) throw "Index out of range in TGlfMultiReader::setActive(const int index)!";
	if(!GLFIsActive[index]){
		GLFIsActive[index] = true;
		activeGLFs.push_back(index);
		pointerToActiveGLFs.push_back(&GLFs[index]);
	}
};

void TGlfMultiReader::_setAllInactive(){
	for(int i=0; i<numGLFs; ++i)
		GLFIsActive[i] = false;
	activeGLFs.clear();
	pointerToActiveGLFs.clear();
};

int TGlfMultiReader::_minChrNumberActiveFiles(){
	int minChr = 9999999;
	for(std::vector<TGlfReader*>::iterator it = pointerToActiveGLFs.begin(); it != pointerToActiveGLFs.end(); ++it){
		if((*it)->chrNumber() < minChr){
			minChr = (*it)->chrNumber();
		}
	}
	return minChr;
};

void TGlfMultiReader::_setCurChrName(){
	//find active file on cur chromosome number and set name
	for(std::vector<TGlfReader*>::iterator it = pointerToActiveGLFs.begin(); it != pointerToActiveGLFs.end(); ++it){
		if((*it)->chrNumber() == _curChrNumber){
			_curChrName = (*it)->chr();
			break;
		}
	}
};

void TGlfMultiReader::_prepareParsing(){
	numActiveFiles = pointerToActiveGLFs.size();
	for(std::vector<TGlfReader*>::iterator it = pointerToActiveGLFs.begin(); it != pointerToActiveGLFs.end(); ++it)
		(*it)->rewind();

	//start at first chromosome, position 0 (one before first position).
	_curChrNumber = 0;
	_curChrLength = -1;
	_position = 1;

	//initialize data
	if(dataInitialized){
		delete[] data;
		delete[] hasData;
	}
	data = new uint8_t*[numActiveFiles];
	hasData = new bool[numActiveFiles];
	dataInitialized = true;
}

void TGlfMultiReader::setActive(const int index){
	if(index >= numGLFs) throw "Index out of range in TGlfMultiReader::setActiveOnly(const int index)!";
	_setAllInactive();
	_setActive(index);
	_prepareParsing();
};

void TGlfMultiReader::setActive(const std::string & name){
	int index = _getGLFIndexFromName(name);
	setActive(index);
};

void TGlfMultiReader::setActive(const int index1, const int index2){
	_setAllInactive();
	_setActive(index1);
	_setActive(index2);
	_prepareParsing();
};

void TGlfMultiReader::setActive(const std::string & name1, const std::string & name2){
	int index1 = _getGLFIndexFromName(name1);
	int index2 = _getGLFIndexFromName(name2);
	setActive(index1, index2);
};

void TGlfMultiReader::setActive(std::vector<int> & indexes){
	_setAllInactive();
	for(std::vector<int>::iterator it=indexes.begin(); it!=indexes.end(); ++it)
		_setActive(*it);
	_prepareParsing();
};

void TGlfMultiReader::setActive(std::vector<std::string> & names){
	_setAllInactive();
	for(std::vector<std::string>::iterator it=names.begin(); it!=names.end(); ++it){
		int index = _getGLFIndexFromName(*it);
		_setActive(index);
	}
	_prepareParsing();
};

void TGlfMultiReader::setAllActive(){
	activeGLFs.clear();
	for(int i=0; i<numGLFs; ++i)
		_setActive(i);
	_prepareParsing();
};

//-------------------------------------
//Looping over active files
//-------------------------------------
bool TGlfMultiReader::moveToNextChromosome(){
	//increment chromosome number
	++_curChrNumber;

	//advance all active files behind in chromosome number
	bool allFilesReachedEnd = true;
	for(std::vector<TGlfReader*>::iterator it = pointerToActiveGLFs.begin(); it != pointerToActiveGLFs.end(); ++it){
		while(!(*it)->eof() && (*it)->chrNumber() < _curChrNumber)
			(*it)->jumpToNextChr();
		if(!(*it)->eof()) allFilesReachedEnd = false;
	}

	//check if we reached end of all files
	if(allFilesReachedEnd) return false;

	//get name and length from first active file
	_position = 1;
	_curChrName = pointerToActiveGLFs[0]->getNameOfParsedChr(_curChrNumber);
	_curChrLength = pointerToActiveGLFs[0]->getLengthOfParsedChr(_curChrNumber);

	//check that all files share the same name and length for this chromosome
	for(std::vector<TGlfReader*>::iterator it = pointerToActiveGLFs.begin(); it != pointerToActiveGLFs.end(); ++it){
		if((*it)->chrNumber() == _curChrNumber){
			if(_curChrName != (*it)->chr())
				throw "Chrosomome names differ between files '" + pointerToActiveGLFs[0]->name() + "' and '" + (*it)->name() + "'!";
			if(_curChrLength != (*it)->chrLength())
				throw "Chrosomome names differ between files '" + pointerToActiveGLFs[0]->name() + "' and '" + (*it)->name() + "'!";
		}
	}

	return true;
};

bool TGlfMultiReader::readNext(){
	//advance position counter
	++_position;
	if(_position > _curChrLength){
		if(!moveToNextChromosome()) return false;
	}

	//advance all files behind position
	_numActiveFilesWithData = 0;
	int i=0;
	for(std::vector<TGlfReader*>::iterator it = pointerToActiveGLFs.begin(); it != pointerToActiveGLFs.end(); ++it, ++i){
		if(!(*it)->eof() && (*it)->chrNumber() == _curChrNumber && (*it)->position < _position)
			(*it)->readNext();
		if(!(*it)->eof() && (*it)->chrNumber() == _curChrNumber && (*it)->position == _position){
			data[i] = (*it)->genotypeQualities;
			hasData[i] = true;
			++_numActiveFilesWithData;
		} else {
			data[i] = genotypeQualitiesMissingData;
			hasData[i] = false;
		}
	}

	//filter sites (i.e., jump to next)

	return true;
};

void TGlfMultiReader::print(){
	std::cout << std::endl << "Multi Reader at position " << _position << " on chromosome '" << _curChrName << std::endl;
	for(int i=0; i<numActiveFiles; ++i){
		std::cout << "File " << i << ":";
		for(int g=0; g<10; ++g) std::cout << "\t" << unsigned(data[i][g]);
		std::cout << std::endl;
	}
};

void TGlfMultiReader::writeSampleNamesOfActiveFiles(gz::ogzstream & out, std::string sep){
	//sample names are file names without glf ending
	if(numActiveFiles > 0){
		for(int i=0; i<numActiveFiles; ++i){
			std::string name = GLFs[i].name();
			out << sep << readBeforeLast(name, ".glf");
		}
	}
};

void TGlfMultiReader::writeVCFHeader(gz::ogzstream & vcf){
	//make sure the header matches the format used in writeSiteToVCF
	vcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	vcf << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n";
	vcf << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
	vcf << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled normalized genotype likelihoods\">\n";

	//also write header with sample names
	vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	writeSampleNamesOfActiveFiles(vcf, "\t");
	vcf << '\n';

};

void TGlfMultiReader::writeSiteToVCF(gz::ogzstream & vcf, int variantQuality, int & refHomIndex, int & hetIndex, int & altHomIndex, TRandomGenerator & randomGenerator){

	PASS ALLELIC CONBINATION INSTEAD!

	//TODO: find way to harmonize code with MLE caller in TSite
	//write position
	vcf << _curChrName << '\t' << _position <<"\t.\t";

	//write major and minor
	vcf << genoMap.baseToChar[genoMap.genotypeToBase[refHomIndex][0]] << '\t' << genoMap.baseToChar[genoMap.genotypeToBase[altHomIndex][0]] << '\t';

	//write quality of variant
	vcf << variantQuality << '\t';

	//write filter, info and format
	vcf << "\t.\t.\tGT:GQ:DP:PL";

	//now write active samples
	for(int i=0; i<numActiveFiles; ++i){
		if(hasData[i]){
			//find min qual
			int minQual = data[i][refHomIndex];
			if(data[i][hetIndex] < minQual) minQual = data[i][hetIndex];
			if(data[i][altHomIndex] < minQual) minQual = data[i][altHomIndex];

			//get all genotypes with minQual (=MLE)
			std::vector<int> mleGenotypes;
			if(data[i][refHomIndex] == minQual) mleGenotypes.push_back(0);
			if(data[i][hetIndex] == minQual) mleGenotypes.push_back(1);
			if(data[i][altHomIndex] == minQual) mleGenotypes.push_back(2);

			//write MLE genoytpe
			int mleGeno = mleGenotypes[randomGenerator.pickOne(mleGenotypes.size())];
			if(mleGeno == 0) vcf << "\t0/0;";
			else if(mleGeno == 1) vcf << "\t0/1;";
			else vcf << "\t1/1;";

			//write genotype quality
			if(mleGenotypes.size() > 1) vcf << "0;";
			else {
				//find second highest quality
				int secondLowestQual = 999;
				if(data[i][refHomIndex] > minQual) secondLowestQual = data[i][refHomIndex];
				if(data[i][hetIndex] > minQual && data[i][hetIndex] < secondLowestQual) secondLowestQual = data[i][refHomIndex];
				if(data[i][altHomIndex] == minQual && data[i][hetIndex] < secondLowestQual) secondLowestQual = data[i][refHomIndex];
				vcf << round(secondLowestQual - minQual) << ";";
			}

			//write depth
			vcf << GLFs[i].depth << ';';

			//write likelihoods
			vcf << data[i][refHomIndex] - minQual << "," << data[i][hetIndex] - minQual << "," << data[i][altHomIndex] - minQual;

		} else {
			vcf << "\t.";
		}
	}

	//end of line
	vcf << '\n';


};
