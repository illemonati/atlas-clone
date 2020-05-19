/*
 * TChromosomes.cpp
 *
 *  Created on: Apr 25, 2019
 *      Author: linkv
 */

#include "TChromosomes.h"
#include "TLog.h"

//---------------------------------------------------------
// TChromosomes
//---------------------------------------------------------
TChromosomes::TChromosomes(BamTools::SamHeader* BamHeader){
	readChromosomes(BamHeader);
};

void TChromosomes::readChromosomes(BamTools::SamHeader* BamHeader){
	//make sure object is empty
	_chromosomes.clear();

	//copy from BamHeader
	uint16_t num = 0;
	for(BamTools::SamSequenceIterator chrIt=BamHeader->Sequences.Begin(); chrIt!=BamHeader->Sequences.End(); ++chrIt, ++num){
		_chromosomes.emplace_back(_chromosomes.size(), chrIt->Name, stringToLong(chrIt->Length));
	}

	_curChr = _chromosomes.end();

	if(_chromosomes.empty()){
		throw "No chromosomes present in BAM header!";
	}
};

TChromosome& TChromosomes::_find(const std::string chrName) const{
	for(auto& c : _chromosomes){
		if(c.name == chrName)
			return c;
	}
	throw "Chromosome '" + chrName + "' not found in BAM header!";
};


void TChromosomes::limitChr(const std::string limitName){
	_find(limitName);

	bool use = true;
	for(auto& it : _chromosomes){
		if(it.name == limitName){
			use = false;
		}
		it.inUse = use;
	}
};

void TChromosomes::useSpecifiedChr(const std::vector<std::string> & chrNames){
	//set all chromosomes to false
	for(auto& c : _chromosomes){
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
	for(auto& it : _chromosomes){
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
	_curChr = _chromosomes.begin();
};

bool TChromosomes::next(){
	++_curChr;
	return _curChr != _chromosomes.end();
};

bool TChromosomes::end() const{
	return _curChr == _chromosomes.end();
};

void TChromosomes::jumpToEnd(){
	_curChr = _chromosomes.end();
};

//getters
uint16_t TChromosomes::size() const{
	return _chromosomes.size();
}

uint32_t TChromosomes::referenceLength() const{
	int chrNum = 0;
	long totLength = 0;
	for(auto& c : _chromosomes){
		if(c.inUse)
			totLength += c.length;
	}
	return totLength;
};

bool TChromosomes::exists(const std::string chrName) const{
	for(auto& c : _chromosomes){
		if(c.name == chrName)
			return true;
	}
	return false;
};

TChromosome& TChromosomes::getChromosome(const std::string chrName){
	return _find(chrName);
};

TChromosome& TChromosomes::curChromosome(){
	return *_curChr;
};


uint16_t TChromosomes::refID(const std::string chrName) const{
	return _find(chrName).refID;
};

uint16_t TChromosomes::curRefID() const{
	return _curChr->refID;
};


uint32_t TChromosomes::length(const uint16_t index) const{
	return _chromosomes[index].length;
};

uint32_t TChromosomes::curLength() const{
	return _curChr->length;
};


std::string TChromosomes::name(const uint16_t index) const{
	return _chromosomes[index].name;
};

std::string TChromosomes::curName() const{
	return _curChr->name;
};


bool TChromosomes::inUse(const uint16_t index) const{
	return _chromosomes[index].inUse;
};

bool TChromosomes::curInUse() const{
	return _curChr->inUse;
};


uint8_t TChromosomes::ploidy(const uint16_t index) const{
	return _chromosomes[index].ploidy;
};

uint8_t TChromosomes::curPloidy() const{
	return _curChr->ploidy;
};

