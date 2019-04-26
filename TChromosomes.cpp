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
	numChromosomes = -1;
	curChrNumber = -1;
};

TChromosomes::TChromosomes(BamTools::SamHeader* BamHeader){
	bamHeader = BamHeader;
	numChromosomes = bamHeader->Sequences.Size();
	int num = 0;
	for(BamTools::SamSequenceIterator chrIt=bamHeader->Sequences.Begin(); chrIt!=bamHeader->Sequences.End(); ++chrIt, ++num){
		names.push_back(chrIt->Name);
		lengths.push_back(stringToLong(chrIt->Length));
		inUse.push_back(true);
		nameMap.emplace(chrIt->Name, num);
	}

	curChrNumber = -1;
	curChrIterator = bamHeader->Sequences.End();
};

void TChromosomes::limitChr(std::string & limitName){
	if(nameMap.find(limitName) == nameMap.end())
		throw "Chromosome limit not found in BAM header!";
	int limitIndex = nameMap.find(limitName)->second;
	for(int i = limitIndex + 1; i<numChromosomes; ++i){
		inUse[i] = false;
	}
}

void TChromosomes::useSpecifiedChr(std::vector<std::string> & chrNames, TLog* logfile){
	logfile->startIndent("Will limit analysis to the following chromosomes:");

	//set all chromosomes to false
	for(int i=0; i<numChromosomes; ++i)
			inUse[i] = false;

	for(std::vector<std::string>::iterator it=chrNames.begin(); it!=chrNames.end(); ++it){
		//find chromosome
		if(nameMap.find(*it) == nameMap.end()){
			throw "Chromosome " + *it + " not found in BAM header!";
		}
		inUse[nameMap.find(*it)->second] = true;
		logfile->list(*it);
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

long TChromosomes::calcReferenceLength(){
	int chrNum = 0;
	long totLength = 0;
	for(BamTools::SamSequenceIterator chrIterator = bamHeader->Sequences.Begin(); chrIterator!=bamHeader->Sequences.End(); ++chrIterator, ++chrNum)
	    if(inUse[chrNum]) totLength += stringToLong(chrIterator->Length);
	return totLength;
};

//getters

int TChromosomes::getNumberOfChr(){
	return numChromosomes;
}

int TChromosomes::getCurIndex(){
	return curChrNumber;
};

long TChromosomes::getCurLength(){
	//TODO: do we need to check range?
	return lengths[curChrNumber];
};

std::string TChromosomes::getCurName(){
	return names[curChrNumber];
};

std::string TChromosomes::getNameAtIndex(int & index){
	return names[index];
}

int TChromosomes::getIndexFromName(const std::string & chrName){
	if(nameMap.find(chrName) != nameMap.end()){
		return nameMap.find(chrName)->second;
	} else {
		throw "Chromosome " + chrName + " was not found in BAM header!";
	}
}

bool TChromosomes::chrInUse(int & index){
	return inUse[index];
}

bool TChromosomes::curChrInUse(){
	return inUse[curChrNumber];
}
