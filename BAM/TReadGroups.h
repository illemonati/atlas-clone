/*
 * TReadGroups.h
 *
 *  Created on: Nov 27, 2015
 *      Author: wegmannd
 */

#ifndef TREADGROUPS_H_
#define TREADGROUPS_H_

#include "stringFunctions.h"
#include "../bamtools/api/BamReader.h"
#include "../bamtools/api/SamSequenceDictionary.h"
#include "TLog.h"
#include "TParameters.h"
#include <vector>
#include <algorithm>

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
//TReadGroups
//---------------------------------------------------------------
struct readGroup{
public:
	std::string name;
	uint16_t id;
	BamTools::SamReadGroup* object;
};

//---------------------------------------------------------------
//TReadGroups
//---------------------------------------------------------------

class TReadGroups{
private:
	readGroup* _groups;
	uint16_t _numGroups;
	bool _initialized;
	bool _limitReadGroups;
	bool* _inUse;

public:
	TReadGroups();
	~TReadGroups();

	void fill(BamTools::SamHeader & bamHeader);
	uint16_t find(const std::string & name) const;
	uint16_t find(BamTools::BamAlignment & alignment) const;
	bool readGroupExists(const std::string & name) const;
	bool readGroupInUse(const uint16_t & readGroupId) const;
	bool readGroupInUse(const std::string name) const;
	bool readGroupInUse(BamTools::BamAlignment & alignment) const;

	const std::string& getName (const uint16_t readGroupId) const;
	uint16_t size() const;
	void filterReadGroups(std::string readGroupList);
	void printReadgroupsInUse(TLog* logfile) const;
	int addTruncatedOrMergedRG(BamTools::SamHeader & bamHeader, std::string oldReadGroupName, std::string newReadGroupName);
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

#endif /* TREADGROUPS_H_ */
