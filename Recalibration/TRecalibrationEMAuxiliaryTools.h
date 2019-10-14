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

//--------------------------------------------------------------------
// TRecalibrationEMQualityPositionMap
// Look-up tables for position and quality. Only indexes will be stored.
//--------------------------------------------------------------------
class TRecalibrationEMQualityPositionMap{
public:
	int maxQual;
	int maxPos;
	double* eta;
	double* etaSquared;
	double* position;
	double* positionSquared;
	bool initialized;

	TRecalibrationEMQualityPositionMap(){
		initialized = false;
		initialize(500, 255); //TODO: think about default!
	};

	~TRecalibrationEMQualityPositionMap(){
		clear();
	};

	void clear(){
		if(initialized){
			delete[] eta;
			delete[] etaSquared;
			delete[] position;
			delete[] positionSquared;
			initialized = false;
		}
	};

	void initialize(int MaxQual, int MaxPos){
		clear();
		maxQual = MaxQual;
		maxPos = MaxPos;
		eta = new double[maxQual+1];
		etaSquared = new double[maxQual+1];
		position = new double[maxPos+1];
		positionSquared = new double[maxPos+1];
		initialized = true;

		//fill qualities. Use TQualityMap for conversion
		TQualityMap qualiMap;
		for(int q=0; q<=maxQual; q++){
			double eps = qualiMap.qualityToError(q);
			if(eps < 0.0000000001) eps = 0.0000000001;
			else if(eps > 0.9999999999) eps = 0.9999999999;

			eta[q] = log(eps / (1.0 - eps));
			etaSquared[q] = eta[q] * eta[q];
		}

		//fill positions
		for(int p = 0; p<=maxPos; p++){
			position[p] = p;
			positionSquared[p] = p * p;
		}
	};
};


//--------------------------------------------------------------------
// TRecalibrationEMReadData
// Per site data storage
//--------------------------------------------------------------------
struct TRecalibrationEMReadData{
	uint8_t quality;
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

	//DEBUG
	int qualContext[100][25];
	int posContext[300][25];
	int posQual[300][100];

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
