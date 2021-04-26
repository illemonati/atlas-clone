/*
 * TReadGroups.h
 *
 *  Created on: Nov 27, 2015
 *      Author: wegmannd
 */

#ifndef TREADGROUPS_H_
#define TREADGROUPS_H_

#include "stringFunctions.h"
#include "TLog.h"
#include "TParameters.h"
#include <vector>
#include <algorithm>
#include <set>
#include "counters.h"

namespace BAM{

//---------------------------------------------------------------
//TReadGroupMaxLength
//---------------------------------------------------------------
struct TReadGroupMaxLength{
public:
	int maxLen;
	uint16_t truncatedOrMergedReadGroupID;
	std::string truncatedOrMergedReadGroup;
	int sequencingType; //0 = single, 1 = mixed, 2 = paired

	TReadGroupMaxLength(int MaxLen, int TruncatedOrMergedReadGroupID, std::string & TruncatedOrMergedReadGroup, int Type){
		maxLen = MaxLen;
		truncatedOrMergedReadGroupID = TruncatedOrMergedReadGroupID;
		truncatedOrMergedReadGroup = TruncatedOrMergedReadGroup;
		sequencingType = Type;
	};
};

//---------------------------------------------------------------
//TReadGroup
//---------------------------------------------------------------
class TReadGroup{
public:
	uint16_t id; //internal ID

	std::string name_ID;
	mutable std::string barcodeSequence_BC;
	mutable std::string sequencingCenter_CN;
	mutable std::string description_DS;
	mutable std::string productionDate_DT;
	mutable std::string flowOrder_FO;
	mutable std::string keySequence_KS;
	mutable std::string library_LB;
	mutable std::string program_PG;
	mutable std::string predictedInsertSize_PI;
	mutable std::string sequencingTechnology_PL;
	mutable std::string platformModel_PM;
	mutable std::string platformUnit_PU;
	mutable std::string sample_SM;

    //flags
    mutable bool inUse; 						//read groups not in use are ignored when reading
    mutable bool writeToHeader;                 //is false if read group is not in use or replaced by new one

    TReadGroup(const uint16_t ID, const std::string Name);
    TReadGroup(const TReadGroup & other) = default;
    TReadGroup* getPointer(){ return this; };
    std::string compileSamHeader() const;

    bool operator<(const TReadGroup & right) const;
    bool operator<(const std::string & right) const;
};

bool operator<(const std::string & left, const TReadGroup & right);

//---------------------------------------------------------------
//TReadGroups
//---------------------------------------------------------------
class TReadGroups{
private:
	std::set<TReadGroup, std::less<>> _readGroups;
	std::vector<const TReadGroup*> _readGroupsById;
	bool _limitReadGroups;

	void _fillLookupFromId();

public:
	TReadGroups();
	~TReadGroups(){};

	void clear();
	const TReadGroup& add(const std::string name);
	const TReadGroup& addAlternativeRG(const std::string Name, const std::string original);
	uint16_t size() const;
	bool empty() const;

	uint16_t getId(const std::string & name) const;
	const std::string& getName (const uint16_t readGroupId) const;
	const TReadGroup& getReadGroup(const std::string & name);
	bool readGroupExists(const std::string & name) const;
	bool readGroupInUse(const uint16_t & readGroupId) const;
	bool readGroupInUse(const std::string name) const;

	//looping over
	std::set<TReadGroup, std::less<>>::iterator begin(){ return _readGroups.begin(); };
	std::set<TReadGroup, std::less<>>::iterator end(){ return _readGroups.end(); };
	std::set<TReadGroup, std::less<>>::iterator cbegin() const{ return _readGroups.cbegin(); };
	std::set<TReadGroup, std::less<>>::iterator cend() const{ return _readGroups.cend(); };

	void filterReadGroups(std::string readGroupList);
	void removeFromHeader(const std::string name);
	void removeFromHeader(const uint16_t readGroupId);
	void printReadgroupsInUse(TLog* logfile) const;
	void fillVectorWithNames(std::vector<std::string> & vec) const;
	std::string compileSamHeader() const;
};

//--------------------------------------------------------------------------------------
//TReadGroupMap
//Maps bam file read group index to internal index, which may differ in case of pooling
//--------------------------------------------------------------------------------------
class TReadGroupMap{
private:
	TReadGroups* _readGroups;

	uint16_t _origNumReadGroups;
	uint16_t _numReadGroups;
	uint16_t* _readGroupMap; //maps read group index to internal index
	std::vector<uint16_t>* _reverseReadGroupMap;

	void _fillWithoutPooling();
	void _fillFromFile(std::string filename, TLog* logfile);
	void _fillReverseMap();
public:
	TReadGroupMap(TReadGroups* ReadGroups);
	TReadGroupMap(TReadGroups* ReadGroups, const std::string filename, TLog* logfile);

	~TReadGroupMap();

	uint16_t getOrigNumReadGroups();
	uint16_t getNumReadGroups();

	uint16_t getIndex(const uint16_t rg);
	uint16_t getIndex(const std::string readGroupName);
	uint16_t operator[](const uint16_t rg);
	void fillNamesOfReadgroups(uint16_t rg, std::vector<std::string> & names);
};

}; //end namespace

#endif /* TREADGROUPS_H_ */
