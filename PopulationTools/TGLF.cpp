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

void TGlfWriter::newChromosome(std::string name, uint32_t length, uint8_t ploidy){
	if(curChr.name != "")
		write(&zero8, sizeof(uint8_t));

	//save cur info
	curChr.update(name, length, ploidy);

	//write new chromosome: length of label, label, length of ref sequence, ploidy
	uint32_t labelLength = curChr.name.size();

	write(&labelLength, sizeof(uint32_t));
	write(curChr.name.c_str(), name.length() * sizeof(char));
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
	position = 0;
	//maxLL = 0;
	depth = 0;
	RMS_mappingQual = 0;
	positionInFile = 0;
	tmpInt8 = 0;
	_lenRead = 0;
	_eof = true;

	genotypeQualitiesMissingData = new uint16_t[curChr.maxNumLikelihoodValues];
	for(int i=0; i<curChr.maxNumLikelihoodValues; ++i)
		genotypeQualitiesMissingData[i] = 0;
};

std::string TGlfReader::getNameOfParsedChr(uint32_t chrNumber){
	if(curChr.number == chrNumber)
		return curChr.name;
	else if(curChr.number > chrNumber)
		return chromosomesAlreadyParsed[chrNumber].name;
	else throw "TGlfReader does not know name of chromosome " + toString(chrNumber) + ": chromosome not yet read!";
};

long TGlfReader::getLengthOfParsedChr(uint32_t chrNumber){
	if(curChr.number == chrNumber)
		return curChr.length;
	else if(curChr.number > chrNumber)
		return chromosomesAlreadyParsed[chrNumber].length;
	else throw "TGlfReader does not know length of chromosome " + toString(chrNumber) + ": chromosome not yet read!";
};

bool TGlfReader::readChr(){
	//store current chromosome name in list of chromosomes parsed
	if(curChr.name != "")
		chromosomesAlreadyParsed.emplace_back(curChr);

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

	//length
	uint32_t length;
	read(&length, sizeof(uint32_t));

	//ploidy
	uint8_t ploidy;
	read(&ploidy, sizeof(uint8_t));

	//update
	curChr.update(name, length, ploidy);

	//set first position = 0
	position = 0;

	return true;
};

bool TGlfReader::chromosomeParsed(std::string & chr){
	for(std::vector< TGlfChromosome >::iterator it=chromosomesAlreadyParsed.begin(); it!=chromosomesAlreadyParsed.end(); ++it){
		if(it->name == chr)
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

	//maxLL and depth (both uint16)
	//read(&maxLL, sizeof(uint16_t));
	read(&depth, sizeof(uint16_t));

	//root mean square of mapping qualities
	read(&RMS_mappingQual, sizeof(uint8_t));
	//RMS_mappingQual = (int) tmpInt8; DO WE NEED THIS??

	//genotype likelihoods
	read(genotypeQualities, curChr.numLikelihoodValues*sizeof(uint16_t));
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
	version.assign(buffer, 4);

	if(version != "GLFA")
		throw "Non-supported GLF version '" + version + "! Currently only version GLFA is supported.";

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

bool TGlfReader::readNextWindow(std::vector<uint16_t*> & genoLikelihoods, std::string chr, long start, long end){
	//Assumes that windows are read in order: no jumping back!
	if(_eof) return false;
	if(chromosomeParsed(chr)) return false;

	//move to correct chromosome
	if(curChr.name != chr){
		while(curChr.name != chr){
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
			memcpy(genoLikelihoods[index], genotypeQualitiesMissingData, curChr.numLikelihoodValues*sizeof(uint16_t));

		memcpy(genoLikelihoods[index], genotypeQualities, curChr.numLikelihoodValues*sizeof(uint16_t));
		++index; ++i;

		//read next record
		if(!readRecordType() || recordType == 0){
			for(; i<end; ++i, ++index)
				memcpy(genoLikelihoods[index], genotypeQualitiesMissingData, curChr.numLikelihoodValues*sizeof(uint16_t));
			readChr();
			break;
		}
		readSNPRecord();
	}

	return true;
};

void TGlfReader::fillGenotypeQualities(uint16_t* destination){
	//assumes pointer points to
	memcpy(destination, genotypeQualities, curChr.numLikelihoodValues*sizeof(uint16_t));
};

void TGlfReader::fillGenotypeLikelihoods(double* destination, TGlfConverter* converter){
	for(int i=0; i<curChr.numLikelihoodValues; ++i){
		destination[i] = converter->toScaledLikelihood(genotypeQualities[i]);
	}
};

//printing
void TGlfReader::printChr(){
	std::cout << "CHROMOSOME: '" << curChr.name << "' of length " << curChr.length << " and ploidy " << (int) curChr.ploidy << "\n";
};

void TGlfReader::printSite(){
	//std::cout << curChr.name << "\t" << position << "\t" << maxLL << "\t" << depth << "\t" << RMS_mappingQual;
	std::cout << curChr.name << "\t" << position << "\t" << depth << "\t" << RMS_mappingQual;
	for(int i=0; i<curChr.numLikelihoodValues; ++i)
		std::cout << "\t" << unsigned(genotypeQualities[i]);
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
	setDepthFilter(params.getParameterIntWithDefault("minDepth", 0), logfile);
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

	minDepth = 0;
};

TGlfMultiReader::~TGlfMultiReader(){
	if(dataInitialized){
		delete[] data;
		delete[] hasData;
	}

	if(readersOpened){
		delete[] GLFIsActive;
		delete[] isHaploid;
		closeGLF();
	}
};

void TGlfMultiReader::_openGLFs(TLog* logfile){
	numGLFs = GLFNames.size();
	GLFIsActive = new bool[numGLFs];
	isHaploid = new bool[numGLFs];

	//open files
	GLFs = new TGlfReader[numGLFs];
	readersOpened = true;
	logfile->startIndent("Opening " + toString(numGLFs) + " GLF files:");
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
	std::string parameter = params.getParameterString("glf");
	//assume that GLF file names are given in a file if string does not contain ".gz"
	if(!stringContains(parameter,".gz")){
		logfile->list("Reading glf input names from file '" + parameter + "'");
		std::ifstream in(parameter.c_str());
		if(!in) throw "Failed to open file '" + parameter + "'!";

		//read file
		std::vector<std::string> vec;
		while(in.good() && !in.eof()){
			std::string line;
			std::getline(in, line);
			fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
			//skip empty lines
			if(vec.size() > 0){
				GLFNames.push_back(vec[0]);
			}
		}
		in.close();
	} else
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

void TGlfMultiReader::setDepthFilter(int MinDepth, TLog* logfile){
	minDepth = MinDepth;
	if(minDepth > 0)
		logfile->list("Will only keep sites with depth >= " + toString(minDepth) + ".");
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
	for(TGlfReader* it : pointerToActiveGLFs){
		if(it->chrNumber() < minChr){
			minChr = it->chrNumber();
		}
	}
	return minChr;
};

void TGlfMultiReader::_setCurChrName(){
	//find active file on cur chromosome number and set name
	for(TGlfReader* it : pointerToActiveGLFs){
		if(it->chrNumber() == _curChrNumber){
			_curChrName = it->chr();
			break;
		}
	}
};

void TGlfMultiReader::_prepareParsing(){
	numActiveFiles = pointerToActiveGLFs.size();
	for(TGlfReader* it : pointerToActiveGLFs)
		it->rewind();

	//start at first chromosome, position 0 (one before first position).
	_curChrNumber = 0;
	_curChrLength = -1;
	_position = 1;

	//initialize data
	if(dataInitialized){
		delete[] data;
		delete[] hasData;
	}
	data = new uint16_t*[numActiveFiles];
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
	for(TGlfReader* it : pointerToActiveGLFs){
		while(!it->eof() && it->chrNumber() < _curChrNumber)
			it->jumpToNextChr();
		if(!it->eof()) allFilesReachedEnd = false;
	}

	//check if we reached end of all files
	if(allFilesReachedEnd) return false;

	//get name and length from first active file not at end
	_position = 1;
	_curChrLength = -1;
	for(TGlfReader* it : pointerToActiveGLFs){
		if(!it->eof() && it->chrNumber() == _curChrNumber){
			_curChrName = it->getNameOfParsedChr(_curChrNumber);
			_curChrLength = it->getLengthOfParsedChr(_curChrNumber);
			break;
		}
	}
	if(_curChrLength < 0) moveToNextChromosome();

	//check that all files share the same name and length for this chromosome
	for(TGlfReader* it : pointerToActiveGLFs){
		if(it->chrNumber() == _curChrNumber){
			if(_curChrName != it->chr())
				throw "Chrosomome names differ between files '" + pointerToActiveGLFs[0]->name() + "' and '" + it->name() + "': '" + _curChrName + "' != '" + it->chr() + "'!";
			if(_curChrLength != it->chrLength())
				throw "Chrosomome names differ between files '" + pointerToActiveGLFs[0]->name() + "' and '" + it->name() + "'!";
		}
	}

	//store whether this chromosome is haploid for each file
	int i = 0;
	for(TGlfReader* it : pointerToActiveGLFs){
		isHaploid[i] = it->chrIsHaploid();
		i++;
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
	for(TGlfReader* it : pointerToActiveGLFs){
		if(!it->eof() && it->chrNumber() == _curChrNumber && it->position < _position)
			it->readNext();
		if(!it->eof() && it->chrNumber() == _curChrNumber && it->position == _position && it->depth > minDepth){
			data[i] = it->genotypeQualities;
			hasData[i] = true;
			++_numActiveFilesWithData;
		} else {
			data[i] = genotypeQualitiesMissingData;
			hasData[i] = false;
		}
		i++;
	}

	//filter sites (i.e., jump to next)
	return true;
};

void TGlfMultiReader::print(){
	std::cout << std::endl << "Multi Reader at position " << _position << " on chromosome '" << _curChrName << std::endl;
	for(int i=0; i<numActiveFiles; ++i){
		std::cout << "File " << i << ":";
		if(isHaploid[i]){
			for(int g=0; g<4; ++g) std::cout << "\t" << unsigned(data[i][g]);
		} else {
			for(int g=0; g<10; ++g) std::cout << "\t" << unsigned(data[i][g]);
		}
		std::cout << std::endl;
	}
};

void TGlfMultiReader::fill(TPopulationLikehoodLocus & storage, const int alleleicCombination){
	storage.resize(numActiveFiles);
	for(int i=0; i<numActiveFiles; ++i){
		storage[i].isHaploid = isHaploid[i];
		storage[i].isMissing = !hasData[i];

		if(isHaploid[i]){
			storage[i].glfLikelihood_0 = data[i][genoMap.alleleicCombinationToBase[alleleicCombination][0]];
			storage[i].glfLikelihood_1 = data[i][genoMap.alleleicCombinationToBase[alleleicCombination][1]];
			storage[i].glfLikelihood_2 = converter.maxValue();
		} else {
			storage[i].glfLikelihood_0 = data[i][genoMap.alleleicCombinationToGenotypes[alleleicCombination][0]];
			storage[i].glfLikelihood_1 = data[i][genoMap.alleleicCombinationToGenotypes[alleleicCombination][1]];
			storage[i].glfLikelihood_2 = data[i][genoMap.alleleicCombinationToGenotypes[alleleicCombination][2]];
		}
	}
};

void TGlfMultiReader::writeSampleNamesOfActiveFiles(gz::ogzstream & out, std::string sep){
	//sample names are file names without glf ending
	if(numActiveFiles > 0){
		for(TGlfReader* it : pointerToActiveGLFs){
			std::string name = it->name();
			out << sep << readBeforeLast(name, ".glf");
		}
	}
};

void TGlfMultiReader::writeVCFHeader(gz::ogzstream & vcf, bool usePhredLikelihoods){
	//make sure the header matches the format used in writeSiteToVCF
	vcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	vcf << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n";
	vcf << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
	if(usePhredLikelihoods)
		vcf << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled normalized genotype likelihoods\">\n";
	else
		vcf << "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Normalized genotype likelihoods\">\n";

	//also write header with sample names
	vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	writeSampleNamesOfActiveFiles(vcf, "\t");
	vcf << '\n';

};

void TGlfMultiReader::writeSiteToVCF(gz::ogzstream & vcf, const int & varianTQuality, const Base major, const Base minor, TRandomGenerator* randomGenerator, const bool & usePhredLikelihoods){
	//Note: we pass hom/het indexes to maintain the major / minor order! Passing the alleleic combination is not enough
	//TODO: find way to harmonize code with Tcaller
	//write position
	vcf << _curChrName << '\t' << _position <<"\t.\t";

	//write major and minor
	vcf << genoMap.baseToChar[major] << '\t' << genoMap.baseToChar[minor] << '\t';

	//write quality of variant
	vcf << varianTQuality;

	//write filter, info and format
	if(usePhredLikelihoods)
		vcf << "\t.\t.\tGT:GQ:DP:PL";
	else
		vcf << "\t.\t.\tGT:GQ:DP:GL";

	//now write active samples
	for(int i=0; i<numActiveFiles; ++i){
		if(isHaploid[i])
			writeHaploidIndividualToVCF(i, vcf, major, minor, randomGenerator, usePhredLikelihoods);
		else
			writeDiploidIndividualToVCF(i, vcf, major, minor, randomGenerator, usePhredLikelihoods);
	}

	//end of line
	vcf << '\n';
};

void TGlfMultiReader::writeDiploidIndividualToVCF(const int ind, gz::ogzstream & vcf, const Base major, const Base minor, TRandomGenerator* randomGenerator, const bool & usePhredLikelihoods){
	if(hasData[ind]){
		//get genotype indeces
		int refHomIndex = genoMap.genotypeMap[major][major];
		int hetIndex = genoMap.genotypeMap[major][minor];
		int altHomIndex = genoMap.genotypeMap[minor][minor];

		//find min qual
		uint16_t minQual = data[ind][refHomIndex];
		if(data[ind][hetIndex] < minQual) minQual = data[ind][hetIndex];
		if(data[ind][altHomIndex] < minQual) minQual = data[ind][altHomIndex];

		//get all genotypes with minQual (=MLE)
		std::vector<int> mleGenotypes;
		if(data[ind][refHomIndex] == minQual) mleGenotypes.push_back(0);
		if(data[ind][hetIndex] == minQual) mleGenotypes.push_back(1);
		if(data[ind][altHomIndex] == minQual) mleGenotypes.push_back(2);

		//write MLE genoytpe
		int mleGeno = mleGenotypes[randomGenerator->pickOne(mleGenotypes.size())];
		if(mleGeno == 0) vcf << "\t0/0:";
		else if(mleGeno == 1) vcf << "\t0/1:";
		else vcf << "\t1/1:";

		//write genotype quality
		if(mleGenotypes.size() > 1) vcf << "0:";
		else {
			//find second highest quality
			int secondLowestQual = converter.maxValue();
			if(data[ind][refHomIndex] > minQual) secondLowestQual = data[ind][refHomIndex];
			if(data[ind][hetIndex] > minQual && data[ind][hetIndex] < secondLowestQual) secondLowestQual = data[ind][refHomIndex];
			if(data[ind][altHomIndex] == minQual && data[ind][hetIndex] < secondLowestQual) secondLowestQual = data[ind][refHomIndex];
			vcf << round(secondLowestQual - minQual) << ":";
		}

		//write depth
		vcf << GLFs[ind].depth << ':';

		//write likelihoods
		if(usePhredLikelihoods){
			vcf << converter.toPhred(data[ind][refHomIndex] - minQual)
				<< "," << converter.toPhred(data[ind][hetIndex] - minQual)
				<< "," << converter.toPhred(data[ind][altHomIndex] - minQual);
		} else {
			//if is to get rid of -0 in output (and having 0 instead). Maybe there is a better way?
			if(data[ind][refHomIndex] == minQual) vcf << "0";
			else vcf << converter.toLog10(data[ind][refHomIndex] - minQual);

			if(data[ind][hetIndex] == minQual) vcf << ",0";
			else vcf << "," << converter.toLog10(data[ind][hetIndex] - minQual);

			if(data[ind][altHomIndex] == minQual) vcf << ",0";
			else vcf << "," << converter.toLog10(data[ind][altHomIndex] - minQual);
		}

	} else {
		vcf << "\t./.:.:.:.";
	}
};

void TGlfMultiReader::writeHaploidIndividualToVCF(int ind, gz::ogzstream & vcf, Base major, Base minor, TRandomGenerator* randomGenerator, const bool & usePhredLikelihoods){
	if(hasData[ind]){
		//find min qual
		int minQual = data[ind][major];
		if(data[ind][minor] < minQual) minQual = data[ind][minor];

		//get all genotypes with minQual (=MLE)
		std::vector<int> mleGenotypes;
		if(data[ind][major] == minQual) mleGenotypes.push_back(major);
		if(data[ind][minor] == minQual) mleGenotypes.push_back(minor);

		//write MLE genoytpe
		int mleGeno = mleGenotypes[randomGenerator->pickOne(mleGenotypes.size())];
		if(mleGeno == major) vcf << "\t0:";
		else vcf << "\t1:";

		//write genotype quality
		if(mleGeno == major) vcf << round(data[ind][minor] - minQual) << ":";
		else  vcf << round(data[ind][major] - minQual) << ":";

		//write depth
		vcf << GLFs[ind].depth << ':';

		//write likelihoods
		if(usePhredLikelihoods){
			vcf << converter.toPhred(data[ind][major] - minQual) << "," << converter.toPhred(data[ind][minor] - minQual);
		} else {
			//if is to get rid of -0 in output (and having 0 instead). Maybe there is a better way?
			if(data[ind][major] == minQual) vcf << "0";
			else vcf << converter.toLog10(data[ind][major] - minQual);

			if(data[ind][minor] == minQual) vcf << ",0";
			else vcf << "," << converter.toLog10(data[ind][minor] - minQual);
		}
	} else {
		vcf << "\t.:.:.:.";
	}
};

