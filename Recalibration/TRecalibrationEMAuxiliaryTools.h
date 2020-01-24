/*
 * TRecalibrationEMAuxiliaryTools.h
 *
 *  Created on: Mar 11, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMAUXILIARYTOOLS_H_
#define TRECALIBRATIONEMAUXILIARYTOOLS_H_

#include "../TQualityMap.h"
#include "../TGenotypeMap.h"
#include "../stringFunctions.h"
#include "../TBase.h"
#include "../TLog.h"
#include "../TReadGroups.h"



//--------------------------------------------------------------
// TRecalibrationEMTransformationMap
//--------------------------------------------------------------
class TRecalibrationEMTransformationMap{
protected:
	uint16_t size;
	uint16_t max;
	double* map;
	bool initialized;

public:
	TRecalibrationEMTransformationMap(){
		initialized = false;
		clear();
	};

	TRecalibrationEMTransformationMap(const uint16_t Max){
		initialized = false;
		initialize(Max);
	};

	~TRecalibrationEMTransformationMap(){
		clear();
	};

	void clear(){
		if(initialized){
			delete[] map;
			initialized = false;
		}
		max = 0;
		size = 0;
		map = nullptr;
	};

	void initialize(const uint16_t Max){
		max = Max;
		size = Max + 1;
		map = new double[size];
		for(int i=0; i<size; i++){
			map[i] = 0.0;
		}
	};

	void set(const uint16_t x, const double value){
		map[x] = value;
	};

	bool checkRange(const uint16_t val){
		if(val <= max) return true;
		else return false;
	};

	double operator[](const int x){
		return map[x];
	};

	double get(const int x){
		return map[x];
	};
};

//--------------------------------------------------------------------
// TRecalibrationEMQualityPositionMap
// Look-up tables for position and quality. Only indexes will be stored.
//--------------------------------------------------------------------
class TRecalibrationEMQualityErrorMap:public TRecalibrationEMTransformationMap{
private:
	void _fill(){
		TQualityMap qualiMap;
		for(int q=0; q<=max; q++){
			double eps = qualiMap.qualityToError(q);

			if(eps < 0.0000000001) eps = 0.0000000001;
			else if(eps > 0.9999999999) eps = 0.9999999999;

			map[q] = log(eps / (1.0 - eps));
		}
	};

public:
	TRecalibrationEMQualityErrorMap(){
		initialized = false;
		initialize(255);
		_fill();
	};
};

//--------------------------------------------------------------------
// TRecalibrationEMReadData
// Per site data storage
//--------------------------------------------------------------------
struct TRecalibrationEMReadData{
	uint8_t quality;
	uint16_t position;
	float D[4];
	uint8_t context;
	uint16_t readGroup;
	bool isSecond;
	uint16_t fragmentLength;

	void setD(Base base, double PMD_CT, double PMD_GA);
};

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
// Object to store for which qualities and positions data is available.
//--------------------------------------------------------------------
class TRecalibrationEMDataTable{
private:
	unsigned int counts;
	bool countsAssembled;

	void assembleCounts();

public:
	int maxQual;
	uint16_t maxPos;
	unsigned int qualities;

	TRecalibrationEMDataTable(const int MaxQual);
	~TRecalibrationEMDataTable();
	void clear();
	void add(TRecalibrationEMReadData & data);
	int size();
	void fillVectorWithUsedQualities(std::vector<uint16_t> & Q);
};

class TRecalibrationEMDataTables{
public:
	int numReadGroups;
	int maxQual;
	TRecalibrationEMDataTable** tables; //tables[readGroup][first/second]
	unsigned int totalCounts;

	TRecalibrationEMDataTables(const int NumReadGroups, const int MaxQual);
	~TRecalibrationEMDataTables();

	void clear();
	void add(TRecalibrationEMReadData & data);
	void assembleCountsPerReadGroup();
	void fillVectorWithUsedQualities(const int readGroupId, const bool isSecondMate, std::vector<uint16_t> & Q);
	TRecalibrationEMDataTable* getTable(const int readGroupId, const bool isSecondMate);
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
