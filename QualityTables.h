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

//---------------------------------------------------------------
//TQualityTable
//---------------------------------------------------------------
class TQualityTable{
private:
	int maxQ;
	int maxQPlusOne;
	long* counts;
	double* freqs;
	bool initialized;

	//tmp
	long sum;

public:
	TQualityTable();
	TQualityTable(int MaxQ);
	TQualityTable(TQualityTable && other);

	TQualityTable& operator=(TQualityTable && other){
		if(this != &other){
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
		}

		return *this;
	};

	~TQualityTable();

	void init(int MaxQ);
	int getMaxQ();
	void add(const int & qual);
	void add(int* qual, int & len);
	void add(TQualityTable & other);
	long at(int qual);
	void calcFrequencies();
	void write(const std::string & filename);
};

//---------------------------------------------------------------
//TQualityTransformTable
//---------------------------------------------------------------
class TQualityTransformTable{
public:
	int maxQInTable;
	int maxQInTablePlusOne;
	double** table; //old qual / new qual
	bool initialized;

	TQualityTransformTable(int maxPhredIntInTable);
	TQualityTransformTable();

	~TQualityTransformTable();

	void initialize(int maxPhredIntInTable);
	void add(const int oldQuality, const int newQuality);
	double size();
	double printTableReturnRSquared(const std::string filename);
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
