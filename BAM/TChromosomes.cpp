/*
 * TChromosomes.cpp
 *
 *  Created on: Apr 25, 2019
 *      Author: linkv
 */

#include "TChromosomes.h"
#include "TLog.h"

TChromosomes::TChromosomes(){};

TChromosomes::TChromosomes(BamTools::SamHeader* BamHeader){
	readChromosomes(BamHeader);
};

void TChromosomes::readChromosomes(BamTools::SamHeader* BamHeader){
	//make sure object is empty
	chromosomes.clear();

	//copy from BamHeader
	uint16_t num = 0;
	for(BamTools::SamSequenceIterator chrIt=BamHeader->Sequences.Begin(); chrIt!=BamHeader->Sequences.End(); ++chrIt, ++num){
		chromosomes.emplace_back(chrIt->Name, stringToLong(chrIt->Length));
	}

	curChr = chromosomes.end();

	if(chromosomes.empty()){
		throw "No chromosomes present in BAM header!";
	}
};

TChromosome& TChromosomes::_find(const std::string chrName){
	for(auto& c : chromosomes){
		if(c.name == chrName)
			return c;
	}
	throw "Chromosome '" + chrName + "' not found in BAM header!";
};


void TChromosomes::limitChr(const std::string limitName){
	_find(limitName);

	bool use = true;
	for(auto& it : chromosomes){
		if(it.name == limitName){
			use = false;
		}
		it.inUse = use;
	}
};

void TChromosomes::useSpecifiedChr(const std::vector<std::string> & chrNames){
	//set all chromosomes to false
	for(auto& c : chromosomes){
		c.inUse = false;
	}

	//set matching to to true
	for(auto& it : chrNames){
		//find chromosome
		auto c = _find(it);
		c.inUse = true;
	}
};

void TChromosomes::writeUsedChromosomes(TLog* logfile){
	for(auto& it : chromosomes){
		if(it.inUse){
			logfile->list(it.name);
		}
	}

	logfile->endIndent();
};



void TChromosomes::specifyPloidy(std::ifstream & ploidyFile, TLog* logfile){
	logfile->startIndent("Setting ploidy for following chromosomes to:");
	while(ploidyFile.good() && !ploidyFile.eof()){
		//read line
		std::string line;
		std::getline(ploidyFile, line);
		std::vector<std::string> vec;
		fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);

		//skip empty lines
		if(vec.size() > 0){
			auto c = _find(vec[0]);
			c.ploidy = stringToInt(vec[1]);

			//report
			if(c.ploidy == 1){
				logfile->list(c.name + ": haploid");
			} else if(c.ploidy == 2){
				logfile->list(c.name + ": diploid");
			} else {
				throw "Unsupported ploidy '" + toString(c.ploidy) + "'! Currently only support haploid (1) and diploid (2).";
			}
		}
	}
	ploidyFile.close();
	logfile->endIndent();
}

void TChromosomes::setToHaploid(std::vector<std::string> chrNames, TLog* logfile){
	logfile->startIndent("Setting the following chromosomes to be haploid:");
	for(auto& it : chrNames){
		auto c = _find(it);
		c.ploidy = 1;
		logfile->list(c.name);
	}
	logfile->endIndent();
}

//move
void TChromosomes::begin(){
	curChr = chromosomes.begin();
};

bool TChromosomes::next(){
	++curChr;
	return curChr != chromosomes.end();
};

bool TChromosomes::end() const{
	return curChr == chromosomes.end();
};

void TChromosomes::jumpToEnd(){
	curChr = chromosomes.end();
};

uint32_t TChromosomes::referenceLength() const{
	int chrNum = 0;
	long totLength = 0;
	for(BamTools::SamSequenceIterator chrIterator = bamHeader->Sequences.Begin(); chrIterator!=bamHeader->Sequences.End(); ++chrIterator, ++chrNum)
	    if(isInUse[chrNum]) totLength += stringToLong(chrIterator->Length);
	return totLength;
};

//getters
uint16_t TChromosomes::size() const{
	return numChromosomes;
}

uint16_t TChromosomes::curIndex() const{
	return curChrNumber;
};

uint32_t TChromosomes::length(const uint16_t index) const{
	return lengths[index];
};

uint32_t TChromosomes::curLength() const{
	return lengths[curChrNumber];
};

std::string TChromosomes::curName() const{
	return names[curChrNumber];
};

std::string TChromosomes::name(const uint16_t index) const{
	return names[index];
};

uint16_t TChromosomes::getIndexFromName(const std::string chrName) const{
	if(nameMap.find(chrName) != nameMap.end()){
		return nameMap.find(chrName)->second;
	} else {
		throw "Chromosome " + chrName + " was not found in BAM header!";
	}
};

bool TChromosomes::exists(const std::string chrName) const{
	if(nameMap.find(chrName) == nameMap.end()){
		return false;
	} else {
		return true;
	}
};

bool TChromosomes::inUse(const uint16_t index) const{
	return isInUse[index];
};

bool TChromosomes::curInUse() const{
	return isInUse[curChrNumber];
};

uint8_t TChromosomes::ploidy(const uint16_t index) const{
	return ploidies[index];
};

uint8_t TChromosomes::curPloidy() const{
	return ploidies[curChrNumber];
};

