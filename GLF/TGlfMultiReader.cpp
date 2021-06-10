/*
 * TGlfMultireader.cpp
 *
 *  Created on: Feb 11, 2020
 *      Author: wegmannd
 */

#include "TGlfMultiReader.h"

namespace GLF{

using coretools::str::toString;

//----------------------------------------------------
//TMultiGLFData
//----------------------------------------------------
TMultiGLFData::TMultiGLFData(const uint32_t & Size){
	samples.resize(Size);
};

void TMultiGLFData::resize(const uint32_t & Size){
	samples.resize(Size);
};

void TMultiGLFData::fill(PopulationTools::TPopulationLikehoodLocus & storage, const BAM::AllelicCombination & alleleicCombination){
	storage.resize(samples.size());
	for(uint32_t i=0; i<samples.size(); ++i){
		storage[i].isHaploid = samples[i].isHaploid;
		storage[i].isMissing = !samples[i].hasData;

		if(samples[i].isHaploid){
			storage[i].glfLikelihood_0 = samples[i].genotypeLikelihoodsGLF[(uint8_t) alleleicCombination.firstAllele()];
			storage[i].glfLikelihood_1 = samples[i].genotypeLikelihoodsGLF[(uint8_t) alleleicCombination.secondAllele()];
			storage[i].glfLikelihood_2 = BAM::HighPrecisionPhredIntErrorRate::max();
		} else {
			storage[i].glfLikelihood_0 = samples[i].genotypeLikelihoodsGLF[(uint8_t) alleleicCombination.homoFirst()];
			storage[i].glfLikelihood_1 = samples[i].genotypeLikelihoodsGLF[(uint8_t) alleleicCombination.het()];
			storage[i].glfLikelihood_2 = samples[i].genotypeLikelihoodsGLF[(uint8_t) alleleicCombination.homoSecond()];
		}
	}
};

uint32_t TMultiGLFData::totalDepth(){
	uint32_t tot = 0;
	for(uint32_t i=0; i<samples.size(); ++i){
		tot += samples[i].depth;
	}
	return tot;
};

//----------------------------------------------------
//TGlfMultiReaderVcf
//----------------------------------------------------
TGlfMultiReaderVcf::TGlfMultiReaderVcf(const std::string filename, const std::string source, std::vector<std::string> & sampleNames, TRandomGenerator* RandomGenerator){
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

void TGlfMultiReaderVcf::_setMajorMinor(const BAM::Base & refAllele, const BAM::Base & altAllele){
	_ref = refAllele; _alt = altAllele;
	_refHom = BAM::Genotype(_ref,_ref);
	_het = BAM::Genotype(_ref,_alt);
	_altHom = BAM::Genotype(_alt,_alt);
};

void TGlfMultiReaderVcf::writeLikelihood(const BAM::HighPrecisionPhredIntErrorRate & likGlf){
	if(_usePhredScaledLikelihoods){
		vcf << (BAM::PhredErrorRate) likGlf;
	} else {
		if(likGlf == 0) vcf << '0';
		else vcf << (coretools::Log10Probability) likGlf;
	}
};

void TGlfMultiReaderVcf::writeSite(const std::string & chrName, const uint32_t & position, const int & varianTQuality, TMultiGLFData & data, const BAM::Base refAllele, const BAM::Base altAllele){
	//Note: we pass hom/het indexes to maintain the major / minor order! Passing the alleleic combination is not enough
	//TODO: find way to harmonize code with TCaller
	_setMajorMinor(refAllele, altAllele);

	//write position
	vcf << chrName << '\t' << position + 1 <<"\t.\t"; //internal position is zero-based!

	//write ref and alt alleles
	vcf << refAllele << '\t' << altAllele;

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
	for(uint32_t i=0; i<data.size(); ++i){
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
		BAM::HighPrecisionPhredIntErrorRate minQual = sample.genotypeLikelihoodsGLF[(uint8_t) _refHom];
		if(sample.genotypeLikelihoodsGLF[(uint8_t) _het] < minQual) minQual = sample.genotypeLikelihoodsGLF[(uint8_t) _het];
		if(sample.genotypeLikelihoodsGLF[(uint8_t) _altHom] < minQual) minQual = sample.genotypeLikelihoodsGLF[(uint8_t) _altHom];

		//get all genotypes with minQual (=MLE)
		std::vector<uint8_t> mleGenotypes;
		if(sample.genotypeLikelihoodsGLF[(uint8_t) _refHom] == minQual) mleGenotypes.push_back(2);
		if(sample.genotypeLikelihoodsGLF[(uint8_t) _het] == minQual) mleGenotypes.push_back(3);
		if(sample.genotypeLikelihoodsGLF[(uint8_t) _altHom] == minQual) mleGenotypes.push_back(4);

		//write MLE genoytpe
		int mleGeno = mleGenotypes[randomGenerator->sample(mleGenotypes.size())];
		vcf << '\t' << genotypeStrings[mleGeno] << ':';

		//write genotype quality
		if(mleGenotypes.size() > 1) vcf << "0:";
		else {
			//find second highest quality
			BAM::HighPrecisionPhredIntErrorRate secondLowestQual = BAM::HighPrecisionPhredIntErrorRate::max();
			if(sample.genotypeLikelihoodsGLF[(uint8_t) _refHom] > minQual){
				secondLowestQual = sample.genotypeLikelihoodsGLF[(uint8_t) _refHom];
			}
			if(sample.genotypeLikelihoodsGLF[(uint8_t) _het] > minQual && sample.genotypeLikelihoodsGLF[(uint8_t) _het] < secondLowestQual){
				secondLowestQual = sample.genotypeLikelihoodsGLF[(uint8_t) _het];
			}
			if(sample.genotypeLikelihoodsGLF[(uint8_t) _altHom] > minQual && sample.genotypeLikelihoodsGLF[(uint8_t) _altHom] < secondLowestQual){
				secondLowestQual = sample.genotypeLikelihoodsGLF[(uint8_t) _altHom];
			}
			vcf << (BAM::PhredIntErrorRate) (secondLowestQual - minQual) << ":";
		}

		//write depth
		vcf << sample.depth << ':';

		//write likelihoods
		writeLikelihood(sample.genotypeLikelihoodsGLF[(uint8_t) _refHom] - minQual);
		vcf << ',';
		writeLikelihood(sample.genotypeLikelihoodsGLF[(uint8_t) _het] - minQual);
		vcf << ',';
		writeLikelihood(sample.genotypeLikelihoodsGLF[(uint8_t) _altHom] - minQual);
	} else {
		vcf << "\t./.:.:.:.";
	}
};

void TGlfMultiReaderVcf::writeHaploidIndividualToVCF(TMultiGLFDataSample & sample){
	if(sample.hasData){
		//find min qual
		BAM::HighPrecisionPhredIntErrorRate minQual = sample.genotypeLikelihoodsGLF[(uint8_t) _ref];
		if(sample.genotypeLikelihoodsGLF[(uint8_t) _alt] < minQual) minQual = sample.genotypeLikelihoodsGLF[(uint8_t) _alt];

		//get all genotypes with minQual (=MLE)
		std::vector<BAM::Base> mleGenotypes;
		if(sample.genotypeLikelihoodsGLF[(uint8_t) _ref] == minQual) mleGenotypes.push_back(_ref);
		if(sample.genotypeLikelihoodsGLF[(uint8_t) _alt] == minQual) mleGenotypes.push_back(_alt);

		//write MLE genoytpe
		BAM::Base mleGeno = mleGenotypes[randomGenerator->sample(mleGenotypes.size())];
		vcf << '\t' << genotypeStrings[mleGeno.get()] << ':';

		//write genotype quality
		if(mleGeno == _ref) vcf << (BAM::PhredIntErrorRate) (sample.genotypeLikelihoodsGLF[(uint8_t) _alt] - minQual) << ":";
		else  vcf << (BAM::PhredIntErrorRate) (sample.genotypeLikelihoodsGLF[(uint8_t) _ref] - minQual) << ":";

		//write depth
		vcf << sample.depth << ':';

		//write likelihoods
		writeLikelihood(sample.genotypeLikelihoodsGLF[(uint8_t) _ref] - minQual);
		vcf << ',';
		writeLikelihood(sample.genotypeLikelihoodsGLF[(uint8_t) _alt] - minQual);
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
	setDepthFilter(params.getParameterWithDefault<int>("minDepth", 0), logfile);
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
		genotypeQualitiesMissingData[i] = BAM::HighPrecisionPhredIntErrorRate::min();

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
	std::string parameter = params.getParameter<std::string>("glf");
	//assume that GLF file names are given in a file if string does not contain ".gz"
	if(!coretools::str::stringContains(parameter,".gz")){
		logfile->list("Reading glf input names from file '" + parameter + "'");
		std::ifstream in(parameter.c_str());
		if(!in) throw "Failed to open file '" + parameter + "'!";

		//read file
		std::vector<std::string> vec;
		while(in.good() && !in.eof()){
			std::string line;
			std::getline(in, line);
			coretools::str::fillContainerFromStringWhiteSpace(line, vec, true);
			//skip empty lines
			if(vec.size() > 0){
				GLFNames.push_back(vec[0]);
			}
		}
		in.close();
	} else {
		params.fillParameterIntoContainer("glf", GLFNames, ',');
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
	fastaBuffer.initialize(FastaFile, 1000000);
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
			for(int g=0; g<4; ++g) std::cout << "\t" << data.samples[i].genotypeLikelihoodsGLF[g];
		} else {
			for(int g=0; g<10; ++g) std::cout << "\t" << data.samples[i].genotypeLikelihoodsGLF[g];
		}
		std::cout << std::endl;
	}
};

void TGlfMultiReader::writeSampleNamesOfActiveFiles(gz::ogzstream & out, std::string sep){
	//sample names are file names without glf ending
	if(numActiveFiles > 0){
		for(TGlfReader* it : pointerToActiveGLFs){
			std::string name = it->name();
			out << sep << coretools::str::readBeforeLast(name, ".glf");
		}
	}
};

void TGlfMultiReader::fillSampleNamesOfActiveFiles(std::vector<std::string> & vec){
	vec.clear();

	//sample names are file names without glf ending
	if(numActiveFiles > 0){
		for(TGlfReader* it : pointerToActiveGLFs){
			vec.emplace_back(coretools::str::readBeforeLast(it->name(), ".glf"));
		}
	}
};

BAM::Base TGlfMultiReader::refBase(){
	if(hasReference){
		return fastaBuffer.refAt(BAM::TGenomePosition(_curRefId, _position));
	} else return BAM::N;
};


}; //end namespace GLF
