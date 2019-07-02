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
	int i;
	long sum;

public:
    TQualityTable();
    TQualityTable(int MaxQ);

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
    void printTable(const std::string filename);
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
    ~TQualityTransformTables();

    void add(const int readGroup, const int oldQuality, const int newQuality);
    void writeTables(std::string outputName);
};

#endif /* QUALITYTABLES_H_ */
