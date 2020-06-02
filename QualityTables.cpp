/*
 * TQualityTables.cpp
 *
 *  Created on: Oct 29, 2019
 *      Author: linkv
 */

/*
 * QualityTables.h
 *
 *  Created on: Nov 30, 2017
 *      Author: phaentu
 */

#include "QualityTables.h"

//---------------------------------------------------------------
//TQualityTable
//---------------------------------------------------------------
TQualityTable::TQualityTable(){
	initialized = false;
	maxPhredInt = -1;
	maxPhredIntPlusOne = -1;
	counts = NULL;
	freqs = NULL;
	sum = 0;
};

TQualityTable::TQualityTable(int MaxPhredInt){
	init(MaxPhredInt);
};

TQualityTable::TQualityTable(TQualityTable && other):maxPhredInt(0),maxPhredIntPlusOne(0),counts(nullptr),freqs(nullptr),initialized(false),sum(0){
	//copy from other
	maxPhredInt = other.maxPhredInt;
	maxPhredIntPlusOne = other.maxPhredIntPlusOne;
	counts = other.counts;
	freqs = other.freqs;
	initialized = other.initialized;
	sum = other.sum;

	//set other to default
	other.maxPhredInt = 0;
	other.maxPhredIntPlusOne = 0;
	other.counts = nullptr;
	other.freqs = nullptr;
	other.initialized = false;
	other.sum = 0;
};


void TQualityTable::init(int MaxPhredInt){
	maxPhredInt = MaxPhredInt; //Note: quality is phred(error) + 33!
	maxPhredIntPlusOne = maxPhredInt + 1;
	counts = new long[maxPhredIntPlusOne];
	for(int i=0; i<maxPhredIntPlusOne; ++i)
		counts[i] = 0;
	freqs = new double[maxPhredIntPlusOne];
	sum = 0;
	initialized = true;
};

TQualityTable::~TQualityTable(){
	if(initialized){
		delete[] counts;
		delete[] freqs;
	}
};


int TQualityTable::getMaxQ(){
	return maxPhredInt;
};

void TQualityTable::add(const int & phredInt){
	if(phredInt < maxPhredIntPlusOne)
		++counts[phredInt];
};

void TQualityTable::add(int* phredInt, int & len){
	for(int i=0; i<len; ++i){
		if(phredInt[i] < maxPhredIntPlusOne)
			++counts[phredInt[i]];
	}
};

void TQualityTable::add(TQualityTable & other){
	int otherMaxQ = other.getMaxQ();
	int m = std::min(maxPhredInt, otherMaxQ) + 1;
	for(int i=33; i<m; ++i)
		counts[i] += other.at(i);
};

long TQualityTable::at(int phredInt){
	return counts[phredInt];
};

void TQualityTable::calcFrequencies(){
	sum = 0;

	for(int i=0; i<maxPhredIntPlusOne; ++i)
		sum += counts[i];

	for(int i=0; i<maxPhredIntPlusOne; ++i)
		freqs[i] = (double) counts[i] / (double) sum;
};

void TQualityTable::write(const std::string & filename){
	//open output file
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";

	//write header
	out << "Quality\tQuality(char)\tCounts\tFrequencies\tCumulativeFrequencies\n";

	//calc frequencies
	calcFrequencies();

	double cumulFreq = 0.0;
	for(int i=0; i<maxPhredIntPlusOne; ++i){
		cumulFreq += freqs[i];
		out << i << "\t" << (char) i << "\t" << counts[i] << "\t" << freqs[i] << "\t" << cumulFreq << "\n";
	}
};



