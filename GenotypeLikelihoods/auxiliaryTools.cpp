/*
 * TRecalibrationEMAuxiliaryTools.cpp
 *
 *  Created on: Mar 11, 2019
 *      Author: phaentu
 */

#include "auxiliaryTools.h"

namespace GenotypeLikelihoods{


//--------------------------------------------------------------------
// TBaseLikelihoods (can also be used as haploid genotype likelihoods)
//--------------------------------------------------------------------
TBaseLikelihoods::TBaseLikelihoods(){
	reset();
};

void TBaseLikelihoods::operator=(const TBaseLikelihoods & other){
	likelihoods[A] = other.at(A);
	likelihoods[C] = other.at(C);
	likelihoods[G] = other.at(G);
	likelihoods[T] = other.at(T);
	likelihoods[N] = other.at(N);
};

void TBaseLikelihoods::reset(){
	likelihoods[A] = 1.0;
	likelihoods[C] = 1.0;
	likelihoods[G] = 1.0;
	likelihoods[T] = 1.0;
	likelihoods[N] = 1.0;
};

//--------------------------------------------------------------------
// TGenotypeLikelihoods
//--------------------------------------------------------------------
TGenotypeLikelihoods::TGenotypeLikelihoods(){
	reset();
};

void TGenotypeLikelihoods::operator=(const TGenotypeLikelihoods & other){
	likelihoods[AA] = other.at(AA);
	likelihoods[AC] = other.at(AC);
	likelihoods[AG] = other.at(AG);
	likelihoods[AT] = other.at(AT);
	likelihoods[CC] = other.at(CC);
	likelihoods[CG] = other.at(CG);
	likelihoods[CT] = other.at(CT);
	likelihoods[GG] = other.at(GG);
	likelihoods[GT] = other.at(GT);
	likelihoods[TT] = other.at(TT);
};

void TGenotypeLikelihoods::reset(){
	set(1.0);
};

void TGenotypeLikelihoods::set(const double val){
	likelihoods[AA] = val;
	likelihoods[AC] = val;
	likelihoods[AG] = val;
	likelihoods[AT] = val;
	likelihoods[CC] = val;
	likelihoods[CG] = val;
	likelihoods[CT] = val;
	likelihoods[GG] = val;
	likelihoods[GT] = val;
	likelihoods[TT] = val;
};



void TGenotypeLikelihoods::fill(const std::vector<TBaseLikelihoods> & bases){
	fill(bases, bases.size());
};

void TGenotypeLikelihoods::fill(const std::vector<TBaseLikelihoods> & bases, const size_t size){
	//allows for vector to be longer than what is to be used
	//do in log if depth is high
	if(bases.size() > 50){
		//initialize
		set(0.0);

		//add to log genotype likelihoods
		for(size_t i=0; i<size; ++i){
			likelihoods[AA] += log(bases[i].at(A));
			likelihoods[AC] += log(0.5*bases[i].at(A) + 0.5*bases[i].at(C));
			likelihoods[AG] += log(0.5*bases[i].at(A) + 0.5*bases[i].at(G));
			likelihoods[AT] += log(0.5*bases[i].at(A) + 0.5*bases[i].at(T));
			likelihoods[CC] += log(bases[i].at(C));
			likelihoods[CG] += log(0.5*bases[i].at(C) + 0.5*bases[i].at(G));
			likelihoods[CT] += log(0.5*bases[i].at(C) + 0.5*bases[i].at(T));
			likelihoods[GG] += log(bases[i].at(G));
			likelihoods[GT] += log(0.5*bases[i].at(G) + 0.5*bases[i].at(T));
			likelihoods[TT] += log(bases[i].at(T));
		}

		//standardize and de-log
		double max = *std::max_element(&likelihoods[AA], &likelihoods[TT]);
		likelihoods[AA] = exp(likelihoods[AA] - max);
		likelihoods[AC] = exp(likelihoods[AC] - max);
		likelihoods[AG] = exp(likelihoods[AG] - max);
		likelihoods[AT] = exp(likelihoods[AT] - max);
		likelihoods[CC] = exp(likelihoods[CC] - max);
		likelihoods[CG] = exp(likelihoods[CG] - max);
		likelihoods[CT] = exp(likelihoods[CT] - max);
		likelihoods[GG] = exp(likelihoods[GG] - max);
		likelihoods[GT] = exp(likelihoods[GT] - max);
		likelihoods[TT] = exp(likelihoods[TT] - max);
	} else { //on natural scale
		//initialize
		set(1.0);

		for(size_t i=0; i<size; ++i){
			likelihoods[AA] *= bases[i].at(A);
			likelihoods[AC] *= 0.5*bases[i].at(A) + 0.5*bases[i].at(C);
			likelihoods[AG] *= 0.5*bases[i].at(A) + 0.5*bases[i].at(G);
			likelihoods[AT] *= 0.5*bases[i].at(A) + 0.5*bases[i].at(T);
			likelihoods[CC] *= bases[i].at(C);
			likelihoods[CG] *= 0.5*bases[i].at(C) + 0.5*bases[i].at(G);
			likelihoods[CT] *= 0.5*bases[i].at(C) + 0.5*bases[i].at(T);
			likelihoods[GG] *= bases[i].at(G);
			likelihoods[GT] *= 0.5*bases[i].at(G) + 0.5*bases[i].at(T);
			likelihoods[TT] *= bases[i].at(T);
		}
	}
};

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
TRecalibrationEMDataTable::TRecalibrationEMDataTable(){
	maxQual = 0;
	maxFragmentLength = 0;
	maxMQ = 0;
	maxPos = 0;
	counts = 0;
	initialized = false;
	countsAssembled = false;
	qualities = nullptr;
	fragmentLengths = nullptr;
	MQ = nullptr;
};

TRecalibrationEMDataTable::TRecalibrationEMDataTable(const int MaxQual, const int MaxFragmentLength, const int maxMQ){
	initialize(MaxQual, MaxFragmentLength, maxMQ);
};

TRecalibrationEMDataTable::~TRecalibrationEMDataTable(){
	delete[] qualities;
	delete[] fragmentLengths;
	delete[] MQ;
};

void TRecalibrationEMDataTable::initialize(const int MaxQual, const int MaxFragmentLength, const int MaxMQ){
	if(initialized){
		delete[] qualities;
		delete[] fragmentLengths;
		delete[] MQ;
	}
	maxQual = MaxQual;
	qualities = new unsigned int[maxQual];

	maxFragmentLength = MaxFragmentLength;
	fragmentLengths = new unsigned int[maxFragmentLength];

	maxMQ = MaxMQ;
	MQ = new unsigned int[maxMQ];

	initialized = true;
	clear();
};

void TRecalibrationEMDataTable::add(TRecalibrationEMReadData & data){
	/*
	if(data.quality > maxQual){
		throw "Can not add data point to TRecalibrationEMDataTable: quality > maxQual!";
	}
	*/
	++qualities[data.qualityAsPhredInt];
	++fragmentLengths[data.fragmentLength];
	++MQ[data.mappingQuality];
	if(maxPos < data.positionFrom5Prime)
		maxPos = data.positionFrom5Prime;
};

void TRecalibrationEMDataTable::assembleCounts(){
	if(!countsAssembled){
		counts = 0;
		for(int i=0; i<maxQual; ++i){
			counts += qualities[i];
		}
		countsAssembled = true;
	}
};

void TRecalibrationEMDataTable::clear(){
	maxPos = 0;
	for(int q=0; q<maxQual; ++q)
		qualities[q] = 0;
	counts = 0;
	countsAssembled = false;
};

size_t TRecalibrationEMDataTable::size(){
	assembleCounts();
	return counts;
}

void TRecalibrationEMDataTable::fillVectorWithUsedQualities(std::vector<uint16_t> & Q){
	Q.clear();
	for(int i=0; i<maxQual; ++i){
		if(qualities[i] > 0){
			Q.push_back(i);
		}
	};
};

void TRecalibrationEMDataTable::fillVectorWithUsedFragmentLengths(std::vector<uint16_t> & lengths){
	lengths.clear();
	for(int i=0; i<maxFragmentLength; ++i){
		if(fragmentLengths[i] > 0){
			lengths.push_back(i);
		}
	};
};

void TRecalibrationEMDataTable::fillVectorWithUsedMQ(std::vector<uint16_t> & MQ){
	MQ.clear();
	for(int i=0; i<maxMQ; ++i){
		if(MQ[i] > 0){
			MQ.push_back(i);
		}
	};
};

//--------------------------------------------------------------------

TRecalibrationEMDataTables::TRecalibrationEMDataTables(const  int NumReadGroups, const int MaxQual, const int MaxFragmentLength, const int MaxMQ){
	numReadGroups = NumReadGroups;
	maxQual = MaxQual;

	tables = new TRecalibrationEMDataTable*[numReadGroups];
	for(int rg = 0; rg<numReadGroups; ++rg){
		tables[rg] = new TRecalibrationEMDataTable[2];
		tables[rg][0].initialize(maxQual, MaxFragmentLength, MaxMQ);
		tables[rg][1].initialize(maxQual, MaxFragmentLength, MaxMQ) ;
	}

	clear();
};

TRecalibrationEMDataTables::~TRecalibrationEMDataTables(){
	for(int rg = 0; rg<numReadGroups; ++rg){
		delete[] tables[rg];
	}

	delete[] tables;
};

void TRecalibrationEMDataTables::clear(){
	for(int rg = 0; rg<numReadGroups; ++rg){
		tables[rg][0].clear();
		tables[rg][1].clear();
	}
	totalCounts = 0;
};

void TRecalibrationEMDataTables::add(TRecalibrationEMReadData & data){
	tables[data.readGroup][(int) data.isSecond].add(data);
};

void TRecalibrationEMDataTables::assembleCountsPerReadGroup(){
	totalCounts = 0;
	for(int rg = 0; rg<numReadGroups; ++rg){
		totalCounts = tables[rg][0].size() + tables[rg][1].size();
	}
};

void TRecalibrationEMDataTables::fillVectorWithUsedQualities(const int readGroupId, const bool isSecondMate, std::vector<uint16_t> & Q){
	tables[readGroupId][(int) isSecondMate].fillVectorWithUsedQualities(Q);
};

TRecalibrationEMDataTable* TRecalibrationEMDataTables::getTable(const int readGroupId, const bool isSecondMate){
	return &tables[readGroupId][(int) isSecondMate];
};

TRecalibrationEMDataTable* TRecalibrationEMDataTables::getTable(const int readGroupId, const int isSecondMate){
	return &tables[readGroupId][isSecondMate];
};
//--------------------------------------------------------------------
// TRecalibrationEMReadGroupIndex
//--------------------------------------------------------------------
TRecalibrationEMReadGroupIndex::TRecalibrationEMReadGroupIndex(){
	readGroups = nullptr;
	readGroupMap = nullptr;
	_numIndices = 0;
	_numReadGroups = 0;
	_numCases = 0;
	_numCasesWithIndex = 0;
	initialized = false;
	readGroupInUse = nullptr;
	readGroupIndex = nullptr;
};

TRecalibrationEMReadGroupIndex::TRecalibrationEMReadGroupIndex(TReadGroups* ReadGroups, TReadGroupMap* ReadGroupMap){
	readGroups = ReadGroups;
	readGroupMap = ReadGroupMap;
	initialized = false;
	_init();
};


TRecalibrationEMReadGroupIndex::~TRecalibrationEMReadGroupIndex(){
	_free();
};

void TRecalibrationEMReadGroupIndex::_init(){
	_numReadGroups = readGroupMap->getNumReadGroups();
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

void TRecalibrationEMReadGroupIndex::initialize(TReadGroups* ReadGroups, TReadGroupMap* ReadGroupMap){
	readGroups = ReadGroups;
	readGroupMap = ReadGroupMap;
	_init();
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

void TRecalibrationEMReadGroupIndex::setAsNotUsed(int readGroup, bool isSecondMate){
	readGroupInUse[readGroup][isSecondMate] = false;
	--_numIndices;
	--_numCasesWithIndex;
	for(int rg = 0; rg<_numReadGroups; ++rg){
		for(int j=0; j<2; ++j){
			if(readGroupIndex[rg][j] > readGroupIndex[readGroup][isSecondMate]){
				--readGroupIndex[rg][j];
			}
		}
	}
	readGroupIndex[readGroup][isSecondMate] = -1;
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

bool TRecalibrationEMReadGroupIndex::hasCasesWithoutIndex(){
	if(_numCasesWithIndex < _numCases)
		return true;
	return false;
};

bool TRecalibrationEMReadGroupIndex::hasCasesWithoutSecondMate(){
	//check if there are read groups without second mate
	for(size_t r=0; r<readGroups->size(); ++r){
		size_t index = readGroupMap->getIndex(r);
		if(!readGroupInUse[index][1])
			return true;
	}
	return false;
};

void TRecalibrationEMReadGroupIndex::reportReadGroupsNotUsed(TLog* logfile){
	for(size_t r=0; r<readGroups->size(); ++r){
		int index = readGroupMap->getIndex(r);
		if(!readGroupInUse[index][0])
			logfile->list(readGroups->getName(r) + " (first mate)");
		if(!readGroupInUse[index][1])
			logfile->list(readGroups->getName(r) + " (second mate)");
	}
};

void TRecalibrationEMReadGroupIndex::reportReadGroupsConsideredSingleEnd(TLog* logfile){
	for(size_t r=0; r<readGroups->size(); ++r){
		int index = readGroupMap->getIndex(r);
		if(!readGroupInUse[index][1])
			logfile->list(readGroups->getName(r));
	}
};

void TRecalibrationEMReadGroupIndex::warningForMissingReadGroups(TLog* logfile){
	for(size_t r=0; r<readGroups->size(); ++r){
		int index = readGroupMap->getIndex(r);
		if(!readGroupInUse[index][0])
			logfile->warning("No recal parameters provided for " + readGroups->getName(r) + " (first mate of paired or single-end). Are you using the correct recal file?");
	}
};


}; //end namespace
