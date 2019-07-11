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

TChromosomes::TChromosomes(SamHeader* BamHeader){
	bamHeader = BamHeader;
    numChromosomes = bamHeader->GetSequences().Size();
	int num = 0;
    SamSequenceIter* iterator = bamHeader->GetSequences().CreateIterator();
    while(iterator->HasNext()){
        SamSequence& sequence = iterator->Next();
        names.push_back(sequence.GetName());
        lengths.push_back(sequence.GetLength());
        isInUse.push_back(true);
        ploidies.push_back(2);
        nameMap.emplace(sequence.GetName(), num);
        ++num;
    }

	curChrNumber = -1;
    curChrIterator = bamHeader->GetSequences().CreateIterator();
};

void TChromosomes::limitChr(std::string & limitName){
	if(nameMap.find(limitName) == nameMap.end())
		throw "Chromosome limit not found in BAM header!";
	int limitIndex = nameMap.find(limitName)->second;
	for(int i = limitIndex + 1; i<numChromosomes; ++i){
		isInUse[i] = false;
	}
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
    curChrIterator = bamHeader->GetSequences().CreateIterator();
};

void TChromosomes::next(){
	++curChrNumber;
    curChrIterator->Next();
};

bool TChromosomes::end(){
    return !curChrIterator->HasNext();
};

void TChromosomes::jumpToEnd(){
    curChrIterator = bamHeader->GetSequences().CreateIteratorEnd();
	curChrNumber = -1;
};

void TChromosomes::jumpToBeginningOfLastChr(){
    if(bamHeader->GetSequences().Size() > 0){
        curChrIterator = bamHeader->GetSequences().CreateIteratorEnd(1);
		curChrNumber = numChromosomes - 1;
	} else
		throw "Found no chromosomes in BAM header!";
};

long TChromosomes::referenceLength(){
	int chrNum = 0;
	long totLength = 0;
    SamSequenceIter* iterator = bamHeader->GetSequences().CreateIterator();
    while(iterator->HasNext()){
        SamSequence& sequence = iterator->Next();
        if(isInUse[chrNum]) totLength += sequence.GetLength();
    }
	return totLength;
};

//getters
int TChromosomes::size(){
	return numChromosomes;
}

int TChromosomes::curIndex(){
	return curChrNumber;
};

long TChromosomes::length(const int index){
	return lengths[index];
};

long TChromosomes::curLength(){
	//TODO: do we need to check range?
	return lengths[curChrNumber];
};

std::string TChromosomes::curName(){
	return names[curChrNumber];
};

std::string TChromosomes::name(const int index){
	return names[index];
}

int TChromosomes::getIndexFromName(const std::string chrName){
	if(nameMap.find(chrName) != nameMap.end()){
		return nameMap.find(chrName)->second;
	} else {
		throw "Chromosome " + chrName + " was not found in BAM header!";
	}
}

bool TChromosomes::inUse(const int index){
	return isInUse[index];
};

bool TChromosomes::curInUse(){
	return isInUse[curChrNumber];
};

int TChromosomes::ploidy(const int index){
	return ploidies[index];
};

int TChromosomes::curPloidy(){
	return ploidies[curChrNumber];
};

