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
	std::string barcodeSequence_BC;
	std::string sequencingCenter_CN;
	std::string description_DS;
	std::string productionDate_DT;
	std::string flowOrder_FO;
	std::string keySequence_KS;
    std::string library_LB;
    std::string program_PG;
    std::string predictedInsertSize_PI;
    std::string sequencingTechnology_PL;
    std::string platformModel_PM;
    std::string platformUnit_PU;
    std::string sample_SM;

    //flags
    bool inUse; 						//read groups not in use are ignored when reading
    bool writeToHeader;                 //is false if read group is not in use or replaced by new one

    TReadGroup(const uint16_t ID, const std::string Name);
    TReadGroup(const TReadGroup & other);
    TReadGroup* getPointer(){ return this; };
    std::string compileSamHeader() const;

    bool operator<(const TReadGroup & right);
    bool operator<(const std::string & left, const TReadGroup & right);
    bool operator<(const TReadGroup & left, const std::string & right);
};

//---------------------------------------------------------------
//TReadGroups
//---------------------------------------------------------------
class TReadGroups{
private:
	std::set<TReadGroup, std::less<>> _readGroups;
	std::vector<TReadGroup*> _readGroupsById;
	bool _limitReadGroups;

	void _fillLookupFromId();

public:
	TReadGroups();
	~TReadGroups(){};

	void clear();
	TReadGroup& add(const std::string name);
	TReadGroup& addAlternativeRG(const std::string Name, const std::string original);
	uint16_t size() const;

	uint16_t getId(const std::string & name) const;
	const std::string& getName (const uint16_t readGroupId) const;
	TReadGroup& getReadGroup(const std::string & name);
	bool readGroupExists(const std::string & name) const;
	bool readGroupInUse(const uint16_t & readGroupId) const;
	bool readGroupInUse(const std::string name) const;

	//looping over
	std::set<TReadGroup, std::less<>>::iterator begin(){ return _readGroups.begin(); };
	std::set<TReadGroup, std::less<>>::iterator end(){ return _readGroups.end(); };

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
