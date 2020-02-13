/*
 * TChromosomes.cpp
 *
 *  Created on: Apr 25, 2019
 *      Author: linkv
 */

#include "TChromosomes.h"
#include "TLog.h"

TChromosomes::TChromosomes(){
	bamHeader = nullptr;
	numChromosomes = 0;
	curChrNumber = 0;
};

TChromosomes::TChromosomes(BamTools::SamHeader* BamHeader){
	bamHeader = BamHeader;
	numChromosomes = bamHeader->Sequences.Size();
	uint16_t num = 0;
	for(BamTools::SamSequenceIterator chrIt=bamHeader->Sequences.Begin(); chrIt!=bamHeader->Sequences.End(); ++chrIt, ++num){
		names.push_back(chrIt->Name);
		lengths.push_back(stringToLong(chrIt->Length));
		isInUse.push_back(true);
		ploidies.push_back(2);
		nameMap.emplace(chrIt->Name, num);
	}

	curChrNumber = -1;
	curChrIterator = bamHeader->Sequences.End();
};

int TChromosomes::limitChr(std::string & limitName){
	if(nameMap.find(limitName) == nameMap.end())
		throw "Chromosome limit not found in BAM header!";
	int limitIndex = nameMap.find(limitName)->second;
	for(int i = limitIndex + 1; i<numChromosomes; ++i){
		isInUse[i] = false;
	}

	//return index of limiting chromosome
	return limitIndex + 1;
}

void TChromosomes::useSpecifiedChr(std::vector<std::string> & chrNames, TLog* logfile){
	logfile->startIndent("Will limit analysis to the following chromosomes:");

	//set all chromosomes to false
	for(int i=0; i<numChromosomes; ++i)
			isInUse[i] = false;

	for(std::vector<std::string>::iterator it=chrNames.begin(); it!=chrNames.end(); ++it){
		//find chromosome
		if(nameMap.find(*it) == nameMap.end()){
			throw "Chromosome " + *it + " not found in BAM header!";
		}
		isInUse[nameMap.find(*it)->second] = true;
		logfile->list(*it);
	}
	logfile->endIndent();
}

void TChromosomes::specifyPloidy(std::ifstream & ploidyFile, TLog* logfile){
	logfile->startIndent("Setting ploidy for following chromosomes to:");
	while(ploidyFile.good() && !ploidyFile.eof()){
		std::string line;
		std::getline(ploidyFile, line);
		std::vector<std::string> vec;
		fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(nameMap.find(vec[0]) == nameMap.end())
				throw "Chromosome " + vec[0] + " not found in BAM header!";
			ploidies[nameMap.find(vec[0])->second] = stringToInt(vec[1]);
			logfile->list(vec[0] + ": " + toString(ploidies[nameMap.find(vec[0])->second]));
		}
	}
	ploidyFile.close();
	logfile->endIndent();
}

void TChromosomes::setToHaploid(std::vector<std::string> chrNames, TLog* logfile){
	logfile->startIndent("Setting the following chromosomes to be haploid:");
	for(std::vector<std::string>::iterator it=chrNames.begin(); it!=chrNames.end(); ++it){
		if(nameMap.find(*it) == nameMap.end()){
			throw "Chromosome " + *it + " not found in BAM header!";
		}
		logfile->list(*it);
		ploidies[nameMap.find(*it)->second] = 1;
	}
	logfile->endIndent();
}

//move

void TChromosomes::begin(){
	curChrNumber = 0;
	curChrIterator = bamHeader->Sequences.Begin();
};

void TChromosomes::next(){
	++curChrNumber;
	++curChrIterator;
};

bool TChromosomes::end(){
	return curChrIterator == bamHeader->Sequences.End();
};

void TChromosomes::jumpToEnd(){
	curChrIterator = bamHeader->Sequences.End();
	curChrNumber = -1;
};

void TChromosomes::jumpToBeginningOfLastChr(){
	if(bamHeader->Sequences.Size() > 0){
		curChrIterator = bamHeader->Sequences.End() - 1;
		curChrNumber = numChromosomes - 1;
	} else
		throw "Found no chromosomes in BAM header!";
};

uint32_t TChromosomes::referenceLength(){
	int chrNum = 0;
	long totLength = 0;
	for(BamTools::SamSequenceIterator chrIterator = bamHeader->Sequences.Begin(); chrIterator!=bamHeader->Sequences.End(); ++chrIterator, ++chrNum)
	    if(isInUse[chrNum]) totLength += stringToLong(chrIterator->Length);
	return totLength;
};

//getters
uint16_t TChromosomes::size(){
	return numChromosomes;
}

uint16_t TChromosomes::curIndex(){
	return curChrNumber;
};

uint32_t TChromosomes::length(const uint16_t index){
	return lengths[index];
};

uint32_t TChromosomes::curLength(){
	//TODO: do we need to check range?
	return lengths[curChrNumber];
};

std::string TChromosomes::curName(){
	return names[curChrNumber];
};

std::string TChromosomes::name(const uint16_t index){
	return names[index];
}

uint16_t TChromosomes::getIndexFromName(const std::string chrName){
	if(nameMap.find(chrName) != nameMap.end()){
		return nameMap.find(chrName)->second;
	} else {
		throw "Chromosome " + chrName + " was not found in BAM header!";
	}
}

bool TChromosomes::inUse(const uint16_t index){
	return isInUse[index];
};

bool TChromosomes::curInUse(){
	return isInUse[curChrNumber];
};

uint8_t TChromosomes::ploidy(const uint16_t index){
	return ploidies[index];
};

uint8_t TChromosomes::curPloidy(){
	return ploidies[curChrNumber];
};

