/*
 * TRecalibrationEMAuxiliaryTools.cpp
 *
 *  Created on: Mar 11, 2019
 *      Author: phaentu
 */

#include "TRecalibrationEMAuxiliaryTools.h"

//--------------------------------------------------------------------
// TRecalibrationEMReadData
//--------------------------------------------------------------------
void TRecalibrationEMReadData::setD(Base base, double PMD_CT, double PMD_GA){
	switch(base){
		case A: D[0] = 0.0; //geno = AA
				D[1] = 1.0; //geno = CC
				D[2] = 1.0 - PMD_GA; //geno = GG
				D[3] = 1.0; //geno = TT
				break;
		case C: D[0] = 1.0; //geno = AA
				D[1] = PMD_CT; //geno = CC
				D[2] = 1.0; //geno = GG
				D[3] = 1.0; //geno = TT
				break;
		case G: D[0] = 1.0; //geno = AA
				D[1] = 1.0; //geno = CC
				D[2] = PMD_GA; //geno = GG
				D[3] = 1.0; //geno = TT
				break;
		case T: D[0] = 1.0; //geno = AA
				D[1] = 1.0 - PMD_CT; //geno = CC
				D[2] = 1.0; //geno = GG
				D[3] = 0.0; //geno = TT
				break;
		case N:
				D[0] = 0.0;
				D[1] = 0.0;
				D[2] = 0.0;
				D[3] = 0.0;
				break;
	}
};

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
//--------------------------------------------------------------------
TRecalibrationEMDataTable::TRecalibrationEMDataTable(int NumReadGroups, int MaxQual){
	numReadGroups = NumReadGroups;
	maxQual = MaxQual;

	qualities = new bool**[numReadGroups];
	maxPos = new unsigned int*[numReadGroups];
	countsPerReadGroup = new unsigned int*[numReadGroups];
	for(int rg = 0; rg<numReadGroups; ++rg){
		qualities[rg] = new bool*[2];
		qualities[rg][0] = new bool[maxQual];
		qualities[rg][1] = new bool[maxQual];

		countsPerReadGroup[rg] = new unsigned int[2];
		maxPos[rg] = new unsigned int[2];
	}

	clear();
};

TRecalibrationEMDataTable::~TRecalibrationEMDataTable(){
	for(int rg = 0; rg<numReadGroups; ++rg){
		delete[] qualities[rg][0];
		delete[] qualities[rg][1];
		delete[] qualities[rg];

		delete[] countsPerReadGroup[rg];
		delete[] maxPos[rg];
	}

	delete[] qualities;
	delete[] maxPos;
	delete[] countsPerReadGroup;
};

void TRecalibrationEMDataTable::clear(){
	for(int rg = 0; rg<numReadGroups; ++rg){
		for(int i=0; i<maxQual; ++i){
			qualities[rg][0][i] = 0;
			qualities[rg][1][i] = 0;
		}
		maxPos[rg][0] = 0;
		maxPos[rg][1] = 0;
		countsPerReadGroup[rg][0] = 0;
		countsPerReadGroup[rg][1] = 0;
	}
	totalCounts = 0;
};

void TRecalibrationEMDataTable::add(TRecalibrationEMReadData & data){
	++qualities[data.readGroup][(int) data.isSecond][data.quality];
	if(maxPos[data.readGroup][data.isSecond] < data.position)
		maxPos[data.readGroup][data.isSecond] = data.position;
};

void TRecalibrationEMDataTable::assembleCountsPerReadGroup(){
	totalCounts = 0;
	for(int rg = 0; rg<numReadGroups; ++rg){
		countsPerReadGroup[rg][0] = 0;
		countsPerReadGroup[rg][1] = 0;
		for(int i=0; i<maxQual; ++i){
			countsPerReadGroup[rg][0] += qualities[rg][0][i];
			countsPerReadGroup[rg][1] += qualities[rg][1][i];
		}
		totalCounts += countsPerReadGroup[rg][0] + countsPerReadGroup[rg][1];
	}
};

//--------------------------------------------------------------------
// TRecalibrationEMReadGroupIndex
//--------------------------------------------------------------------
TRecalibrationEMReadGroupIndex::TRecalibrationEMReadGroupIndex(){
	_numIndices = 0;
	_numReadGroups = 0;
	_numCases = 0;
	_numCasesWithIndex = 0;
	initialized = false;
	readGroupInUse = NULL;
	readGroupIndex = NULL;
};

TRecalibrationEMReadGroupIndex::~TRecalibrationEMReadGroupIndex(){
	_free();
}

void TRecalibrationEMReadGroupIndex::_init(int NumReadGroups){
	_numReadGroups = NumReadGroups;
	_numCases = 2 * _numReadGroups;
	_numCasesWithIndex = 0;
	_numIndices = 0;
	readGroupInUse = new bool*[_numReadGroups];
	readGroupIndex = new int*[_numReadGroups];

	for(int rg = 0; rg<_numReadGroups; ++rg){
		readGroupInUse[rg] = new bool[2];
		readGroupInUse[rg][0] = false;
		readGroupInUse[rg][1] = false;

		readGroupIndex[rg] = new int[2];
		readGroupIndex[rg][0] = -1;
		readGroupIndex[rg][1] = -1;
	}

	initialized = true;
};

void TRecalibrationEMReadGroupIndex::initialize(int NumReadGroups){
	_init(NumReadGroups);
};

void TRecalibrationEMReadGroupIndex::_free(){
	if(initialized){
		for(int rg = 0; rg<_numReadGroups; ++rg){
			delete[] readGroupInUse[rg];
			delete[] readGroupIndex[rg];
		}

		delete[] readGroupInUse;
		delete[] readGroupIndex;
	}
	initialized = false;
};

void TRecalibrationEMReadGroupIndex::setAllAsUsed(){
	for(int rg = 0; rg<_numReadGroups; ++rg){
		readGroupInUse[rg][0] = true;
		readGroupInUse[rg][1] = true;
		readGroupIndex[rg][0] = 2*rg;
		readGroupIndex[rg][1] = 2*rg + 1;
	}
	_numIndices = 2*_numReadGroups;
	_numCasesWithIndex = 2*_numReadGroups;
};

void TRecalibrationEMReadGroupIndex::setAllToSingleIndex(){
	for(int rg = 0; rg<_numReadGroups; ++rg){
		readGroupInUse[rg][0] = true;
		readGroupInUse[rg][1] = true;
		readGroupIndex[rg][0] = 0;
		readGroupIndex[rg][1] = 0;
	}
	_numIndices = 1;
	_numCasesWithIndex = 2*_numReadGroups;
};

int TRecalibrationEMReadGroupIndex::setAsUsed(int readGroup, bool isSecondMate){
	readGroupInUse[readGroup][isSecondMate] = true;
	readGroupIndex[readGroup][isSecondMate] = _numIndices;
	++_numIndices;
	++_numCasesWithIndex;
	return readGroupIndex[readGroup][isSecondMate];
};

bool TRecalibrationEMReadGroupIndex::nextNotInUse(std::pair<int, bool> & pair){
	//check if there are read groups not in use. If so, return true and fill pair. Otherwise, return false
	for(int rg = 0; rg<_numReadGroups; ++rg){
		for(int j=0; j<2; ++j){
			if(!readGroupInUse[rg][j]){
				pair.first = rg;
				pair.second = j;
				return true;
			}
		}
	}
	return false;
};

void TRecalibrationEMReadGroupIndex::reportReadGroupsNotUsed(TLog* logfile, TReadGroups & readGroups){
	for(int rg = 0; rg<_numReadGroups; ++rg){
		if(!readGroupInUse[rg][0])
			logfile->list(readGroups.getName(rg) + " (first mate)");
		if(!readGroupInUse[rg][1])
			logfile->list(readGroups.getName(rg) + " (second mate)");
	}
};

