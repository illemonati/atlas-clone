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


//---------------------------------------------------------------
//TQualityTransformTable
//---------------------------------------------------------------

TQualityTransformTable::TQualityTransformTable(int MaxPhredInt){
	initialized = false;
	initialize(MaxPhredInt);
};

TQualityTransformTable::TQualityTransformTable(){
	initialized = false;
	maxPhredInt = 0;
	maxPhredIntPlusOne = 0;
	table = NULL;
};

TQualityTransformTable::~TQualityTransformTable(){
	for(int i=0; i<maxPhredIntPlusOne; ++i){
		delete[] table[i];
	}
	delete[] table;
};

void TQualityTransformTable::initialize(int MaxPhredInt){
	if(initialized == true)
		throw "Quality table already initialized!";

	maxPhredInt = MaxPhredInt;
	maxPhredIntPlusOne = maxPhredInt + 1;

	//initialize table
	table = new double*[maxPhredIntPlusOne];
	for(int i=0; i<maxPhredIntPlusOne; ++i){
		table[i] = new double[maxPhredIntPlusOne];
		for(int j=0; j<maxPhredIntPlusOne; ++j){
			table[i][j] = 0;
		}
	}
	initialized = true;
};

void TQualityTransformTable::add(const int oldPhredInt, const int newPhredInt){
	if(oldPhredInt < maxPhredIntPlusOne && newPhredInt < maxPhredIntPlusOne){
		table[oldPhredInt][newPhredInt] += 1.0;
	}
};

double TQualityTransformTable::size(){
	double size = 0;
	for(int i=0; i<maxPhredIntPlusOne; ++i){
		for(int j=0; j<maxPhredIntPlusOne; ++j){
			size += table[i][j];
		}
	}
	return size;
};

double TQualityTransformTable::writeTableAndReturnRSquared(const std::string filename){
	//get total
	double totSum = size();

	//calculate mean oldQ
	double meanOldQ = 0.0;
	for(int i=0; i<maxPhredIntPlusOne; ++i){
		int colSum = 0;
		for(int j=0; j<maxPhredIntPlusOne; ++j){
			colSum += table[i][j];
		}
		meanOldQ += (double) i * (double) colSum;
	}
	meanOldQ /= totSum;

	//open file
	TOutputFilePlain out(filename);

	//print header
	std::vector<std::string> header = {"oldQ/newQ"};
	for(int i=0; i<maxPhredIntPlusOne; ++i)
		header.push_back(toString(i));

	//prepare calculating R squared
	double SSres = 0.0;
	double SStot = 0.0;

	//print rows
	for(int i=0; i<maxPhredIntPlusOne; ++i){
		out << i;
		for(int j=33; j<maxPhredIntPlusOne; ++j){
			double weight = table[i][j] / totSum;
			SSres += weight * ((double) i - (double) j) * ((double) i - (double) j);
			SStot += weight * ((double) j - meanOldQ) * ((double) j - meanOldQ);
		}
		out.endLine();
	}

	return((SStot - SSres) / SStot);
};


//---------------------------------------------------------------
//TQualityTransformTables
//---------------------------------------------------------------
TQualityTransformTables::TQualityTransformTables(TReadGroups & ReadGroups, int MaxPhredInt){
	readGroups = &ReadGroups;

	combinedTable.initialize(MaxPhredInt);
	perReadGroupTables = new TQualityTransformTable[readGroups->size()];
	for(unsigned int i=0; i<readGroups->size(); i++)
		perReadGroupTables[i].initialize(MaxPhredInt);
};

void TQualityTransformTables::add(const int readGroup, const int oldPhredInt, const int newPhredInt){
	perReadGroupTables[readGroup].add(oldPhredInt, newPhredInt);
	combinedTable.add(oldPhredInt, newPhredInt);
};

void TQualityTransformTables::writeTables(std::string outputName, TLog* logfile){
	//print tables for read groups
	for(unsigned int i=0; i<readGroups->size(); ++i){
		double RSquared = perReadGroupTables[i].writeTableAndReturnRSquared(outputName + "_" + readGroups->getName(i) + "_qualityTransformation.txt");
		logfile->conclude("R squared for read group " + readGroups->getName(i) + " is " + toString(RSquared));
	}

	double RSquared = combinedTable.writeTableAndReturnRSquared(outputName + "_total_qualityTransformation.txt");
	logfile->conclude("R squared for total data is " + toString(RSquared));

};


