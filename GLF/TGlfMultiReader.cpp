/*
 * TGlfMultireader.cpp
 *
 *  Created on: Feb 11, 2020
 *      Author: wegmannd
 */

#include "TGlfMultiReader.h"


//----------------------------------------------------
//TMultiGLFData
//----------------------------------------------------
void TMultiGLFData::_free(){
	if(size > 0){
		delete[] samples;
	}
};

TMultiGLFData::TMultiGLFData(){
	size = 0;
	samples = nullptr;
};

TMultiGLFData::TMultiGLFData(int Size){
	size = 0;
	resize(Size);
};

void TMultiGLFData::resize(uint32_t Size){
	if(size != Size){
		_free();
		size = Size;
		if(size > 0){
			samples = new TMultiGLFDataSample[size];
		}
	}
};

void TMultiGLFData::fill(PopulationTools::TPopulationLikehoodLocus & storage, const int alleleicCombination){
	storage.resize(size);
	for(uint32_t i=0; i<size; ++i){
		storage[i].isHaploid = samples[i].isHaploid;
		storage[i].isMissing = !samples[i].hasData;

		if(samples[i].isHaploid){
			storage[i].glfLikelihood_0 = samples[i].genotypeLikelihoodsGLF[genoMap.alleleicCombinationToBase[alleleicCombination][0]];
			storage[i].glfLikelihood_1 = samples[i].genotypeLikelihoodsGLF[genoMap.alleleicCombinationToBase[alleleicCombination][1]];
			storage[i].glfLikelihood_2 = converter.maxValue();
		} else {
			storage[i].glfLikelihood_0 = samples[i].genotypeLikelihoodsGLF[genoMap.alleleicCombinationToGenotypes[alleleicCombination][0]];
			storage[i].glfLikelihood_1 = samples[i].genotypeLikelihoodsGLF[genoMap.alleleicCombinationToGenotypes[alleleicCombination][1]];
			storage[i].glfLikelihood_2 = samples[i].genotypeLikelihoodsGLF[genoMap.alleleicCombinationToGenotypes[alleleicCombination][2]];
		}
	}
};

uint32_t TMultiGLFData::totalDepth(){
	uint32_t tot = 0;
	for(uint32_t i=0; i<size; ++i){
		tot += samples[i].depth;
	}
	return tot;
};

//----------------------------------------------------
//TGlfMultiReaderVcf
//----------------------------------------------------
TGlfMultiReaderVcf::TGlfMultiReaderVcf(const std::string filename, const std::string source, std::vector<std::string> & sampleNames, TRandomGenerator* RandomGenerator){
	ref = N; alt = N;
	ref_char = 'N'; alt_char = 'N';
	refHomIndex = 0; hetIndex = 0; altHomIndex = 0;
	_usePhredScaledLikelihoods = false;
	randomGenerator = RandomGenerator;

	//open vcf file
	_openVCF(filename, source, sampleNames);
};

void TGlfMultiReaderVcf::_openVCF(const std::string & filename, const std::string & source, std::vector<std::string> & sampleNames){
	_closeVCF();

	//open vcf file
	vcf.open(filename.c_str());
	vcfOpened = true;

	//write info
	//TODO: create VCF class to harmonize code across different uses. Also include code in Tiger and other
	vcf << "##fileformat=VCFv4.2\n";
	vcf << "##source=" << source << "\n";

	//make sure the header matches the format used in writeSiteToVCF
	vcf << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
	vcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	vcf << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n";
	vcf << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
	if(_usePhredScaledLikelihoods)
		vcf << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled normalized genotype likelihoods\">\n";
	else
		vcf << "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Normalized genotype likelihoods\">\n";

	//also write header with sample names
	vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	for(std::string& name : sampleNames)
		vcf << '\t' << name;
	vcf << '\n';
};

void TGlfMultiReaderVcf::_closeVCF(){
	vcf.close();
	vcfOpened = false;
};

void TGlfMultiReaderVcf::_setMajorMinor(const Base & refAllele, const Base & altAllele){
	ref = refAllele; alt = altAllele;
	ref_char = genoMap.baseToChar[ref];
	alt_char = genoMap.baseToChar[alt];
	refHomIndex = genoMap.genotypeMap[ref][ref];
	hetIndex = genoMap.genotypeMap[ref][alt];
	altHomIndex = genoMap.genotypeMap[alt][alt];
};

void TGlfMultiReaderVcf::writeLikelihood(const uint16_t & likGlf){
	if(_usePhredScaledLikelihoods){
		vcf << (int) converter.toPhred(likGlf);
	} else {
		if(likGlf == 0) vcf << '0';
		else vcf << converter.toLog10(likGlf);
	}
};

void TGlfMultiReaderVcf::writeSite(const std::string & chrName, const uint32_t & position, const int & varianTQuality, TMultiGLFData & data, const Base refAllele, const Base altAllele){
	//Note: we pass hom/het indexes to maintain the major / minor order! Passing the alleleic combination is not enough
	//TODO: find way to harmonize code with TCaller
	_setMajorMinor(refAllele, altAllele);

	//write position
	vcf << chrName << '\t' << position + 1 <<"\t.\t"; //internal position is zero-based!

	//write ref and alt alleles
	vcf << ref_char << '\t' << alt_char;

	//write quality of variant
	vcf << '\t' << varianTQuality;

	//write info field: total depth
	vcf << "\t.\tDP=" << data.totalDepth();

	//write filter, info and format
	if(_usePhredScaledLikelihoods)
		vcf << "\tGT:GQ:DP:PL";
	else
		vcf << "\tGT:GQ:DP:GL";

	//now write active samples
	for(uint32_t i=0; i<data.size; ++i){
		if(data.samples[i].isHaploid)
			writeHaploidIndividualToVCF(data.samples[i]);
		else
			writeDiploidIndividualToVCF(data.samples[i]);
	}

	//end of line
	vcf << '\n';
};

void TGlfMultiReaderVcf::writeDiploidIndividualToVCF(TMultiGLFDataSample & sample){
	if(sample.hasData){
		//find min qual
		uint16_t minQual = sample.genotypeLikelihoodsGLF[refHomIndex];
		if(sample.genotypeLikelihoodsGLF[hetIndex] < minQual) minQual = sample.genotypeLikelihoodsGLF[hetIndex];
		if(sample.genotypeLikelihoodsGLF[altHomIndex] < minQual) minQual = sample.genotypeLikelihoodsGLF[altHomIndex];

		//get all genotypes with minQual (=MLE)
		std::vector<uint8_t> mleGenotypes;
		if(sample.genotypeLikelihoodsGLF[refHomIndex] == minQual) mleGenotypes.push_back(2);
		if(sample.genotypeLikelihoodsGLF[hetIndex] == minQual) mleGenotypes.push_back(3);
		if(sample.genotypeLikelihoodsGLF[altHomIndex] == minQual) mleGenotypes.push_back(4);

		//write MLE genoytpe
		int mleGeno = mleGenotypes[randomGenerator->sample(mleGenotypes.size())];
		vcf << '\t' << genotypeStrings[mleGeno] << ':';

		//write genotype quality
		if(mleGenotypes.size() > 1) vcf << "0:";
		else {
			//find second highest quality
			int secondLowestQual = converter.maxValue();
			if(sample.genotypeLikelihoodsGLF[refHomIndex] > minQual){
				secondLowestQual = sample.genotypeLikelihoodsGLF[refHomIndex];
			}
			if(sample.genotypeLikelihoodsGLF[hetIndex] > minQual && sample.genotypeLikelihoodsGLF[hetIndex] < secondLowestQual){
				secondLowestQual = sample.genotypeLikelihoodsGLF[hetIndex];
			}
			if(sample.genotypeLikelihoodsGLF[altHomIndex] > minQual && sample.genotypeLikelihoodsGLF[altHomIndex] < secondLowestQual){
				secondLowestQual = sample.genotypeLikelihoodsGLF[altHomIndex];
			}
			vcf << (int) converter.toPhred(secondLowestQual - minQual) << ":";
		}

		//write depth
		vcf << sample.depth << ':';

		//write likelihoods
		writeLikelihood(sample.genotypeLikelihoodsGLF[refHomIndex] - minQual);
		vcf << ',';
		writeLikelihood(sample.genotypeLikelihoodsGLF[hetIndex] - minQual);
		vcf << ',';
		writeLikelihood(sample.genotypeLikelihoodsGLF[altHomIndex] - minQual);
	} else {
		vcf << "\t./.:.:.:.";
	}
};

void TGlfMultiReaderVcf::writeHaploidIndividualToVCF(TMultiGLFDataSample & sample){
	if(sample.hasData){
		//find min qual
		int minQual = sample.genotypeLikelihoodsGLF[ref];
		if(sample.genotypeLikelihoodsGLF[alt] < minQual) minQual = sample.genotypeLikelihoodsGLF[alt];

		//get all genotypes with minQual (=MLE)
		std::vector<int> mleGenotypes;
		if(sample.genotypeLikelihoodsGLF[ref] == minQual) mleGenotypes.push_back(0);
		if(sample.genotypeLikelihoodsGLF[alt] == minQual) mleGenotypes.push_back(1);

		//write MLE genoytpe
		int mleGeno = mleGenotypes[randomGenerator->pickOne(mleGenotypes.size())];
		vcf << '\t' << genotypeStrings[mleGeno] << ':';

		//write genotype quality
		if(mleGeno == ref) vcf << (int) converter.toPhred(sample.genotypeLikelihoodsGLF[alt] - minQual) << ":";
		else  vcf << (int) converter.toPhred(sample.genotypeLikelihoodsGLF[ref] - minQual) << ":";

		//write depth
		vcf << sample.depth << ':';

		//write likelihoods
		writeLikelihood(sample.genotypeLikelihoodsGLF[ref] - minQual);
		vcf << ',';
		writeLikelihood(sample.genotypeLikelihoodsGLF[alt] - minQual);
	} else {
		vcf << "\t.:.:.:.";
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
	GLFs = nullptr;
	GLFIsActive = nullptr;
	numActiveFiles = 0;
	_onlyJumpToPositionsWithData = false;

	_position = 0;
	_nextPosition = 0;
	_numActiveFilesWithData = 0;

	for(int i=0; i<10; ++i)
		genotypeQualitiesMissingData[i] = 0;

	minDepth = 0;
	hasReference = false;
};

TGlfMultiReader::~TGlfMultiReader(){
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
			fillVectorFromStringWhiteSpace(line, vec, true);
			//skip empty lines
			if(vec.size() > 0){
				GLFNames.push_back(vec[0]);
			}
		}
		in.close();
	} else {
		params.fillParameterIntoVector("glf", GLFNames, ',');
	}
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

void TGlfMultiReader::addReference(const std::string FastaFile){
	fastaBuffer.initialize(FastaFile, &genoMap, 1000000);
	hasReference = true;
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

void TGlfMultiReader::_prepareParsing(){
	//prepare data
	numActiveFiles = pointerToActiveGLFs.size();
	data.resize(numActiveFiles);

	for(TGlfReader* it : pointerToActiveGLFs)
		it->rewind();

	//read first SNP record in all active files
	for(TGlfReader* it : pointerToActiveGLFs){
		it->readNext();
	}

	//where to start?
	if(_onlyJumpToPositionsWithData){
		_jumpToNextPositionWithData();
	} else {
		_curRefId = 0;
		_position = 0;
		_nextPosition = 0;

		//fill chromosome info
		_updateChromosomeInfo();
	}
};

bool TGlfMultiReader::_jumpToNextPositionWithData(){
	//find min(refId) and min(pos)
	bool allFilesReachedEnd = true;
	bool first = true;
	_curRefId = pointerToActiveGLFs[0]->refId();
	_position = pointerToActiveGLFs[0]->position();

	for(TGlfReader* it : pointerToActiveGLFs){
		if(!it->eof()){
			allFilesReachedEnd = false;

			if(first){
				_curRefId = it->refId();
				_position = it->position();
				first = false;
			} else {
				if(it->refId() < _curRefId){
					_curRefId = it->refId();
					_position = it->position();
				} else if(it->refId() == _curRefId && it->position() < _position){
					_position = it->position();
				}
			}
		}
	}

	//did we reach end?
	if(allFilesReachedEnd) return false;

	//fill chromosome info
	_updateChromosomeInfo();

	//make sure first record is used
	_nextPosition = _position;
	return true;
};

void TGlfMultiReader::_updateChromosomeInfo(){
	//update curChr
	TGlfChromosome* chr = pointerToActiveGLFs[0]->pointerToChr(_curRefId);
	if(chr == nullptr){
		throw "Chromosome with refId " + toString(_curRefId) + " not present in GLF file '" + pointerToActiveGLFs[0]->name() + "!";
	}
	_curChr.update(*chr);

	//check that all files share the same chromosome info
	chr = nullptr;
	int i=0;
	for(TGlfReader* it : pointerToActiveGLFs){
		if(it->fillPointerToChr(_curRefId, chr)){
			if(chr->name != _curChr.name)
				throw "Chrosomome names differ between files '" + pointerToActiveGLFs[0]->name() + "' and '" + it->name() + "': '" + _curChr.name + "' != '" + chr->name + "'!";
			if(chr->length != _curChr.length)
				throw "Chrosomome lengths differ between files '" + pointerToActiveGLFs[0]->name() + "' and '" + it->name() + "': '" + toString(_curChr.length) + "' != '" + toString(chr->length) + "'!";

			//store ploidy
			if(chr->ploidy == 1) data.samples[i].isHaploid = true;
			else if(chr->ploidy == 2) data.samples[i].isHaploid = false;
			else throw "Ploidy " + toString(chr->ploidy) + " is currently not supported!";
		} else {
			throw "Chromosome with refId " + toString(_curRefId) + " not present in any GLF file provided! Limit to only sites with data?";
		}
		++i;
	}
};

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
	//increment chromosome ref_char id
	++_curRefId;

	//advance all active files behind
	bool allFilesReachedEnd = true;
	for(TGlfReader* it : pointerToActiveGLFs){
		while(!it->eof() && it->refId() < _curRefId)
			it->jumpToNextChr();
		if(!it->eof()) allFilesReachedEnd = false;
	}

	//check if we reached end of all files
	if(allFilesReachedEnd) return false;

	//get name and length from first active file not at end
	if(_onlyJumpToPositionsWithData){
		_jumpToNextPositionWithData();
	} else {
		_position = 0;
		_nextPosition = 0;
	}

	//get pointer to chromosome
	_updateChromosomeInfo();

	return true;
};

bool TGlfMultiReader::readNext(){
	//advance to next position
	if(_nextPosition > _curChr.length){
		if(!moveToNextChromosome()) return false;
	}

	//advance all files behind next position
	_numActiveFilesWithData = 0;
	bool allFilesReachedEnd = true;
	int i=0;
	for(TGlfReader* it : pointerToActiveGLFs){
		while(!it->eof() && it->refId() == _curRefId && it->position() < _nextPosition){
			it->readNext();
		}

		if(!it->eof() && it->position() == _nextPosition && it->refId() == _curRefId){
			data.samples[i].genotypeLikelihoodsGLF = it->pointerToGenotypeLikelihoodsGLF();
			data.samples[i].hasData = true;
			data.samples[i].depth = it->depth();
			++_numActiveFilesWithData;
		} else {
			data.samples[i].genotypeLikelihoodsGLF = genotypeQualitiesMissingData;
			data.samples[i].hasData = false;
			data.samples[i].depth = 0;
		}

		if(!it->eof()) allFilesReachedEnd = false;

		++i;
	}

	//check if we reached end of all files
	if(allFilesReachedEnd) return false;

	//jump?
	if(_onlyJumpToPositionsWithData && _numActiveFilesWithData == 0){
		if(_jumpToNextPositionWithData())
			return readNext();
		else return false;
	}

	//update position
	_position = _nextPosition;
	++_nextPosition;

	return true;
};

void TGlfMultiReader::print(){
	std::cout << std::endl << "Multi Reader at position " << _position << " on chromosome '" << _curChr.name << std::endl;
	for(int i=0; i<numActiveFiles; ++i){
		std::cout << "File " << i << ":";
		if(data.samples[i].isHaploid){
			for(int g=0; g<4; ++g) std::cout << "\t" << unsigned(data.samples[i].genotypeLikelihoodsGLF[g]);
		} else {
			for(int g=0; g<10; ++g) std::cout << "\t" << unsigned(data.samples[i].genotypeLikelihoodsGLF[g]);
		}
		std::cout << std::endl;
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

void TGlfMultiReader::fillSampleNamesOfActiveFiles(std::vector<std::string> & vec){
	vec.clear();

	//sample names are file names without glf ending
	if(numActiveFiles > 0){
		for(TGlfReader* it : pointerToActiveGLFs){
			vec.emplace_back(readBeforeLast(it->name(), ".glf"));
		}
	}
};

Base TGlfMultiReader::refBase(){
	if(hasReference){
		return fastaBuffer.refAt(BAM::TGenomePosition(_curRefId, _position));
	} else return N;
};



