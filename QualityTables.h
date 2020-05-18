/*
 * QualityTables.h
 *
 *  Created on: Nov 30, 2017
 *      Author: phaentu
 */

#ifndef QUALITYTABLES_H_
#define QUALITYTABLES_H_

#include <fstream>
#include "TReadGroups.h"
#include "TFile.h"

//---------------------------------------------------------------
//TQualityTable
//---------------------------------------------------------------
class TQualityTable{
private:
	int maxPhredInt;
	int maxPhredIntPlusOne;
	long* counts;
	double* freqs;
	bool initialized;

	//tmp
	long sum;

public:
	TQualityTable();
	TQualityTable(int MaxPhredInt);
	TQualityTable(TQualityTable && other);

	TQualityTable& operator=(TQualityTable && other){
		if(this != &other){
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
		}

		return *this;
	};

	~TQualityTable();

	void init(int MaxPhredInt);
	int getMaxQ();
	void add(const int & phredInt);
	void add(int* phredInt, int & len);
	void add(TQualityTable & other);
	long at(int phredInt);
	void calcFrequencies();
	void write(const std::string & filename);
};

//---------------------------------------------------------------
//TQualityTransformTable
//---------------------------------------------------------------
class TQualityTransformTable{
public:
	int maxPhredInt;
	int maxPhredIntPlusOne;
	double** table; //old phredInt / new phredInt
	bool initialized;

	TQualityTransformTable(int maxPhredIntInTable);
	TQualityTransformTable();

	~TQualityTransformTable();

	void initialize(int maxPhredIntInTable);
	void add(const int oldPhredInt, const int newPhredInt);
	double size();
	double writeTableAndReturnRSquared(const std::string filename);
};

//---------------------------------------------------------------
//TQualityTransformTables
//---------------------------------------------------------------
class TQualityTransformTables{
private:
	TReadGroups* readGroups;
	TQualityTransformTable* perReadGroupTables;
	TQualityTransformTable combinedTable;

public:
	TQualityTransformTables(TReadGroups & ReadGroups, int MaxQ);

	~TQualityTransformTables(){
		delete[] perReadGroupTables;
	};

	void add(const int readGroup, const int oldQuality, const int newQuality);
	void writeTables(std::string outputName, TLog* logfile);
};

#endif /* QUALITYTABLES_H_ */
