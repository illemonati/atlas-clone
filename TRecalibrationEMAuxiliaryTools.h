/*
 * TRecalibrationEMAuxiliaryTools.h
 *
 *  Created on: Mar 11, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMAUXILIARYTOOLS_H_
#define TRECALIBRATIONEMAUXILIARYTOOLS_H_

#include "TQualityMap.h"
#include "TGenotypeMap.h"
#include "stringFunctions.h"
#include "TBase.h"
#include "TLog.h"
#include "TReadGroups.h"

//--------------------------------------------------------------------
// TRecalibrationEMQualityPositionMap
// Look-up tables for position and quality. Only indexes will be stored.
//--------------------------------------------------------------------
class TRecalibrationEMQualityPositionMap{
public:
	int maxQualPhredInt;
	int maxPos;
	double* eta;
	double* etaSquared;
	double* position;
	double* positionSquared;
	bool initialized;

	TRecalibrationEMQualityPositionMap();
	~TRecalibrationEMQualityPositionMap();

	void clear();
	void initialize(int MaxQualPhredInt, int MaxPos);
};

//--------------------------------------------------------------------
// TRecalibrationEMReadData
// Per site data storage
//--------------------------------------------------------------------
struct TRecalibrationEMReadData{
	uint8_t qualityPhredInt;
	uint8_t position;
	float D[4];
	uint8_t context;
	uint16_t readGroup;
	bool isSecond;

	void setD(Base base, double PMD_CT, double PMD_GA);
};

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
// Object to store for which qualities and positions data is available.
//--------------------------------------------------------------------
class TRecalibrationEMDataTable{
public:
	int numReadGroups;
	int maxQual;
	int	*** qualities; //qualities[readGroup][first/second][quality]
	unsigned int** maxPos; //maxPos[readGroup][first/second]
	unsigned int** countsPerReadGroup;
	unsigned int totalCounts;

	TRecalibrationEMDataTable(int NumReadGroups, int MaxQual);
	~TRecalibrationEMDataTable();

	void clear();
	void add(TRecalibrationEMReadData & data);
	void assembleCountsPerReadGroup();
	void fillVectorWithUsedQualities(int readGroupId, bool isSecondMate, std::vector<int> & Q);
};


//--------------------------------------------------------------------
// TRecalibrationEMReadGroupIndex
// Object to map read group ID and mate to an internal index used to store recal parameters
//--------------------------------------------------------------------
class TRecalibrationEMReadGroupIndex{
private:
	bool initialized;
	int** readGroupIndex;
	bool** readGroupInUse;

	void _init(int NumReadGroups);
	void _free();

	int _numReadGroups;
	int _numCases;  //cases = 2 * numReadGroups as each read group has two mates
	int _numCasesWithIndex;
	int _numIndices;

public:
	TRecalibrationEMReadGroupIndex();
	~TRecalibrationEMReadGroupIndex();

	void initialize(int NumReadGroups);
	int numReadGroups(){ return _numReadGroups; };

	void setAllAsUsed();
	void setAllToSingleIndex();
	int setAsUsed(int readGroup, bool isSecondMate);
	void setAsNotUsed(int readGroup, bool isSecondMate);

	bool inUse(const int readGroup, const bool isSecondMate){
		return readGroupInUse[readGroup][isSecondMate];
	};

	int index(const int readGroup, const bool isSecondMate){
		return readGroupIndex[readGroup][isSecondMate];
	};

	int index(const TRecalibrationEMReadData & data){
		return readGroupIndex[data.readGroup][data.isSecond];
	};

	bool nextNotInUse(std::pair<int, bool> & pair);

	bool hasCasesWithoutIndex(){
		if(_numCasesWithIndex < _numCases)
			return true;
		return false;
	};

	void reportReadGroupsNotUsed(TLog* logfile, TReadGroups & readGroups, TReadGroupMap & ReadGroupMap);
	void reportReadGroupsNotUsed(TLog* logfile, TReadGroups & readGroups);
	void reportReadGroupsConsideredSingleEnd(TLog* logfile, TReadGroups & readGroups, TReadGroupMap & ReadGroupMap);
	void reportReadGroupsConsideredSingleEnd(TLog* logfile, TReadGroups & readGroups);
	void warningForMissingReadGroups(TLog* logfile, TReadGroups & readGroups, TReadGroupMap & ReadGroupMap);
	void warningForMissingReadGroups(TLog* logfile, TReadGroups & readGroups);
};




#endif /* TRECALIBRATIONEMAUXILIARYTOOLS_H_ */
