/*
 * TRecalibrationEMAuxiliaryTools.h
 *
 *  Created on: Mar 11, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMAUXILIARYTOOLS_H_
#define TRECALIBRATIONEMAUXILIARYTOOLS_H_

#include <TBase.h>
#include <TGenotypeData.h>
#include <cstddef>
#include "TQualityMap.h"
#include "TGenotypeMap.h"
#include "stringFunctions.h"
#include "TLog.h"
#include "TReadGroups.h"
#include <algorithm>
#include "TFile.h"

namespace GenotypeLikelihoods{

//--------------------------------------------------------------------
// TRecalibrationEMReadData
// Per site data storage
//--------------------------------------------------------------------
struct TRecalibrationEMReadData{
	uint8_t qualityAsPhredInt;
	uint16_t positionFrom5Prime;
	float D[4];
	uint8_t context;
	uint16_t readGroup;
	bool isSecond;
	uint16_t fragmentLength;
	uint8_t mappingQuality;

	void setD(Base base, double PMD_CT, double PMD_GA);
};

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
// Object to store for which qualities and positions data is available.
//--------------------------------------------------------------------
class TRecalibrationEMDataTable{
private:
	size_t counts;
	bool countsAssembled;
	bool initialized;

	void assembleCounts();

public:
	int maxQual, maxFragmentLength, maxMQ;
	uint16_t maxPos;
	unsigned int* qualities;
	unsigned int* fragmentLengths;
	unsigned int* MQ;

	TRecalibrationEMDataTable();
	TRecalibrationEMDataTable(const int MaxQual, const int MaxFragmentLength, const int MQ);
	~TRecalibrationEMDataTable();

	void initialize(const int MaxQual, const int MaxFragmentLength, const int MQ);
	void clear();
	void add(TRecalibrationEMReadData & data);
	size_t size();
	void fillVectorWithUsedQualities(std::vector<uint16_t> & Q);
	void fillVectorWithUsedFragmentLengths(std::vector<uint16_t> & lengths);
	void fillVectorWithUsedMQ(std::vector<uint16_t> & MQ);
};

class TRecalibrationEMDataTables{
public:
	int numReadGroups;
	int maxQual;
	TRecalibrationEMDataTable** tables; //tables[readGroup][first/second]
	unsigned int totalCounts;

	TRecalibrationEMDataTables(const int NumReadGroups, const int MaxQual, const int MaxFragmentLength, const int MQ);
	~TRecalibrationEMDataTables();

	void clear();
	void add(TRecalibrationEMReadData & data);
	void assembleCountsPerReadGroup();
	void fillVectorWithUsedQualities(const int readGroupId, const bool isSecondMate, std::vector<uint16_t> & Q);
	TRecalibrationEMDataTable* getTable(const int readGroupId, const bool isSecondMate);
	TRecalibrationEMDataTable* getTable(const int readGroupId, const int isSecondMate);
};


//--------------------------------------------------------------------------------------
// TRecalibrationEMReadGroupIndex
// Object to map read group ID and mate to an internal index used to store recal models
//--------------------------------------------------------------------------------------
class TRecalibrationEMReadGroupIndex{
private:
	BAM::TReadGroups* readGroups;
	BAM::TReadGroupMap* readGroupMap;
	bool initialized;
	int** readGroupIndex;
	bool** readGroupInUse;

	void _init();
	void _free();

	int _numReadGroups;
	int _numCases;  //cases = 2 * numReadGroups as each read group has two mates
	int _numCasesWithIndex;
	int _numIndices;

public:
	TRecalibrationEMReadGroupIndex();
	TRecalibrationEMReadGroupIndex(BAM::TReadGroups* ReadGroups, BAM::TReadGroupMap* ReadGroupMap);
	~TRecalibrationEMReadGroupIndex();

	void initialize(BAM::TReadGroups* ReadGroups, BAM::TReadGroupMap* ReadGroupMap);
	int numReadGroups(){ return _numReadGroups; };

	void setAllAsUsed();
	void setAllToSingleIndex();
	int setAsUsed(int readGroup, bool isSecondMate);
	void setAsNotUsed(int readGroup, bool isSecondMate);

	bool inUse(const int readGroup, const bool isSecondMate) const{
		return readGroupInUse[readGroup][isSecondMate];
	};

	int index(const int readGroup, const bool isSecondMate) const{
		return readGroupIndex[readGroup][isSecondMate];
	};

	int index(const TBase & base) const{
		return readGroupIndex[base.readGroupID][base.isSecondMate()];
	};

	int index(const TRecalibrationEMReadData & data) const{
		return readGroupIndex[data.readGroup][data.isSecond];
	};

	//std::string name(int readGroup, bool isSecondMate);
	bool nextNotInUse(std::pair<int, bool> & pair) const;
	bool hasCasesWithoutIndex() const;
	bool hasCasesWithoutSecondMate() const;
	void reportReadGroupsNotUsed(TLog* logfile) const;
	void reportReadGroupsConsideredSingleEnd(TLog* logfile) const;
	void warningForMissingReadGroups(TLog* logfile) const;
};

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


class TRecalibrationEMQualityTransformationMap:public TRecalibrationEMTransformationMap{
public:
	TRecalibrationEMQualityTransformationMap(){
		TQualityMap qualiMap;

		initialize(qualiMap.minPhredInt);

		//now set
		for(uint16_t q=0; q<qualiMap.minPhredInt; ++q){
			//convert phred int quality to error
			double eps = qualiMap.phredIntToError(q);
			if(eps < 0.0000000001) eps = 0.0000000001;
			else if(eps > 0.9999999999) eps = 0.9999999999;

			//then transform into logit space
			map[q] = log(eps / (1.0 - eps));
		}
	}
};

//--------------------------------------------------------------------
// Derivatives
//--------------------------------------------------------------------
struct TRecalibrationEMFirstDerivative{
	uint16_t index;
	double derivative;
};

struct TRecalibrationEMSecondDerivative{
	uint16_t index1;
	uint16_t index2;
	double derivative;
};

typedef std::vector<TRecalibrationEMFirstDerivative>::iterator TRecalibrationEMFirstDerivativesIterator;

class TRecalibrationEMFirstDerivatives{
private:
	std::vector<TRecalibrationEMFirstDerivative> _derivatives;
	TRecalibrationEMFirstDerivativesIterator _cur;

public:
	TRecalibrationEMFirstDerivatives(){};
	TRecalibrationEMFirstDerivatives(size_t Size){
		resize(Size);
	};

	void resize(size_t Size){
		_derivatives.resize(Size);
	};

	size_t size(){ return _derivatives.size(); };
	void  restart(){ _cur = _derivatives.begin(); };

	void add(const uint16_t parameterIndex, const double derivative){
		_cur->index = parameterIndex;
		_cur->derivative = derivative;
		++_cur;
	};

	TRecalibrationEMFirstDerivativesIterator begin(){
		return _derivatives.begin();
	};

	TRecalibrationEMFirstDerivativesIterator end(){
		return _derivatives.end();
	};
};

typedef std::vector<TRecalibrationEMSecondDerivative>::iterator TRecalibrationEMSecondDerivativesIterator;

class TRecalibrationEMSecondDerivatives{
private:
	std::vector<TRecalibrationEMSecondDerivative> _derivatives;
	TRecalibrationEMSecondDerivativesIterator _cur;

public:
	TRecalibrationEMSecondDerivatives(){};
	TRecalibrationEMSecondDerivatives(size_t Size){
		resize(Size);
	};

	void resize(size_t Size){
		_derivatives.resize(Size);
	};

	size_t size(){ return _derivatives.size(); };
	void  restart(){ _cur = _derivatives.begin(); };

	void add(const uint16_t parameterIndex1, const uint16_t parameterIndex2, const double derivative){
		_cur->index1 = parameterIndex1;
		_cur->index2 = parameterIndex2;
		_cur->derivative = derivative;
		++_cur;
	};

	TRecalibrationEMSecondDerivativesIterator begin(){
		return _derivatives.begin();
	};

	TRecalibrationEMSecondDerivativesIterator end(){
		return _derivatives.end();
	};
};


}; //end namespace


#endif /* TRECALIBRATIONEMAUXILIARYTOOLS_H_ */
