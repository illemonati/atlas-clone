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
	maxQ = -1;
	maxQPlusOne = -1;
	counts = NULL;
	freqs = NULL;
	sum = 0;
};

TQualityTable::TQualityTable(int MaxQ){
	init(MaxQ);
};

TQualityTable::TQualityTable(TQualityTable && other):maxQ(0),maxQPlusOne(0),counts(nullptr),freqs(nullptr),initialized(false),sum(0){
	//copy from other
	maxQ = other.maxQ;
	maxQPlusOne = other.maxQPlusOne;
	counts = other.counts;
	freqs = other.freqs;
	initialized = other.initialized;
	sum = other.sum;

	//set other to default
	other.maxQ = 0;
	other.maxQPlusOne = 0;
	other.counts = nullptr;
	other.freqs = nullptr;
	other.initialized = false;
	other.sum = 0;
};


void TQualityTable::init(int MaxQ){
	maxQ = MaxQ; //Note: quality is phred(error) + 33!
	maxQPlusOne = maxQ + 1;
	counts = new long[maxQPlusOne];
	for(int i=33; i<maxQPlusOne; ++i)
		counts[i] = 0;
	freqs = new double[maxQPlusOne];
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
	return maxQ;
};

void TQualityTable::add(const int & qual){
	if(qual < maxQPlusOne)
		++counts[qual];
};

void TQualityTable::add(int* qual, int & len){
	for(int i=0; i<len; ++i){
		if(qual[i] < maxQPlusOne)
			++counts[qual[i]];
	}
};

void TQualityTable::add(TQualityTable & other){
	int otherMaxQ = other.getMaxQ();
	int m = std::min(maxQ, otherMaxQ) + 1;
	for(int i=33; i<m; ++i)
		counts[i] += other.at(i);
};

long TQualityTable::at(int qual){
	return counts[qual];
};

void TQualityTable::calcFrequencies(){
	sum = 0;

	for(int i=33; i<maxQPlusOne; ++i)
		sum += counts[i];

	for(int i=33; i<maxQPlusOne; ++i)
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
	for(int i=33; i<maxQPlusOne; ++i){
		cumulFreq += freqs[i];
		out << i-33 << "\t" << (char) i << "\t" << counts[i] << "\t" << freqs[i] << "\t" << cumulFreq << "\n";
	}
};


//---------------------------------------------------------------
//TQualityTransformTable
//---------------------------------------------------------------

TQualityTransformTable::TQualityTransformTable(int maxPhredIntInTable){
	initialized = false;
	initialize(maxPhredIntInTable);
};

TQualityTransformTable::TQualityTransformTable(){
	initialized = false;
	maxQInTable = 0;
	maxQInTablePlusOne = 0;
	table = NULL;
};

TQualityTransformTable::~TQualityTransformTable(){
	for(int i=0; i<maxQInTablePlusOne; ++i){
		delete[] table[i];
	}
	delete[] table;
};

void TQualityTransformTable::initialize(int maxPhredIntInTable){
	if(initialized == true)
		throw "Quality table already initialized!";

	maxQInTable = maxPhredIntInTable + 33;
	maxQInTablePlusOne = maxQInTable + 1;
	table = new double*[maxQInTablePlusOne];
	for(int i=0; i<maxQInTablePlusOne; ++i){
		table[i] = new double[maxQInTablePlusOne];
		for(int j=0; j<maxQInTablePlusOne; ++j){
			table[i][j] = 0;
		}
	}
	initialized = true;
};

void TQualityTransformTable::add(const int oldQuality, const int newQuality){
	if(oldQuality < maxQInTable && newQuality < maxQInTable){
		table[oldQuality][newQuality] += 1.0;
	}
};

double TQualityTransformTable::size(){
	double size = 0;
	for(int i=33; i<maxQInTablePlusOne; ++i){
		for(int j=33; j<maxQInTablePlusOne; ++j){
			size += table[i][j];
		}
	}
	return size;
};

double TQualityTransformTable::printTableReturnRSquared(const std::string filename){
	//get total
	double sum = size();

	//calculate mean oldQ
	double meanOldQ = 0.0;
	for(int i=33; i<maxQInTablePlusOne; ++i){
		double weight_i = 0.0;
		for(int j=33; j<maxQInTablePlusOne; ++j){
			double weight = table[i][j] / sum;
			weight_i += weight;
		}
		meanOldQ += (double) i * weight_i;
	}

	//open file
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";

	//print header
	out << "oldQ/newQ";
	for(int i=33; i<maxQInTablePlusOne; ++i)
		out << "\t" << i-33;
	out << "\n";

	//prepare calculating R squared
	double SSres = 0.0;
	double SStot = 0.0;

	//print rows
	for(int i=33; i<maxQInTablePlusOne; ++i){
		out << i-33;
		double meanOldQ = 0.0;
		for(int j=33; j<maxQInTablePlusOne; ++j){
			double weight = table[i][j] / sum;
			out << "\t" << weight;
			SSres += weight * ((double) i - (double) j) * ((double) i - (double) j);
			SStot += weight * ((double) j - meanOldQ) * ((double) j - meanOldQ);
		}
		out << "\n";
	}

	//close file
	out.close();

	return((SStot - SSres) / SStot);
};


//---------------------------------------------------------------
//TQualityTransformTables
//---------------------------------------------------------------
TQualityTransformTables::TQualityTransformTables(TReadGroups & ReadGroups, int MaxQ){
	readGroups = &ReadGroups;

	combinedTable.initialize(MaxQ);
	perReadGroupTables = new TQualityTransformTable[readGroups->size()];
	for(unsigned int i=0; i<readGroups->size(); i++)
		perReadGroupTables[i].initialize(MaxQ);
};

void TQualityTransformTables::add(const int readGroup, const int oldQuality, const int newQuality){
	perReadGroupTables[readGroup].add(oldQuality, newQuality);
	combinedTable.add(oldQuality, newQuality);
};

void TQualityTransformTables::writeTables(std::string outputName){
	//print tables for read groups
	for(unsigned int i=0; i<readGroups->size(); ++i)
		perReadGroupTables[i].printTableReturnRSquared(outputName + "_" + readGroups->getName(i) + "_qualityTransformation.txt");

	combinedTable.printTableReturnRSquared(outputName + "_total_qualityTransformation.txt");
};


