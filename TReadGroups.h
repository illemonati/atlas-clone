/*
 * TReadGroups.h
 *
 *  Created on: Nov 27, 2015
 *      Author: wegmannd
 */

#ifndef TREADGROUPS_H_
#define TREADGROUPS_H_

#include "stringFunctions.h"
#include "bamtools/api/BamReader.h"
#include "bamtools/api/SamSequenceDictionary.h"
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
	size_t _numGroups;
	bool _initialized;
	bool _limitReadGroups;
	bool* _inUse;

public:
	TReadGroups();
	~TReadGroups();

	void fill(BamTools::SamHeader & bamHeader);
	int find(const std::string & name);
	int find(BamTools::BamAlignment & alignment);
	bool readGroupExists(const std::string & name);
	bool readGroupInUse(const int & readGroupId);
	bool readGroupInUse(const size_t & readGroupId);
	bool readGroupInUse(const std::string name);
	bool readGroupInUse(BamTools::BamAlignment & alignment);

	std::string getName(int readGroupId);
	size_t size();
	void filterReadGroups(std::string readGroupList);
	void printReadgroupsInUse(TLog* logfile);
	int addTruncatedOrMergedRG(BamTools::SamHeader & bamHeader, std::string oldReadGroupName, std::string newReadGroupName);
};

//--------------------------------------------------------------------------------------
//TReadGroupMap
//Maps bam file read group index to internal index, which may differ in case of pooling
//--------------------------------------------------------------------------------------
class TReadGroupMap{
private:

	TReadGroups* _readGroups;

	int _origNumReadGroups;
	int _numReadGroups;
	int* _readGroupMap; //maps read group index to internal index
	std::vector<int>* _reverseReadGroupMap;

	void _fillWithoutPooling();
	void _fillFromFile(std::string filename, TLog* logfile);
	void _fillReverseMap();
public:
	TReadGroupMap(TReadGroups* ReadGroups);
	TReadGroupMap(TReadGroups* ReadGroups, const std::string filename, TLog* logfile);

	~TReadGroupMap();

	int getOrigNumReadGroups();
	int getNumReadGroups();

	int operator[](int rg);
	int getIndex(int rg);
	int getIndex(const std::string readGroupName);
	void fillNamesOfReadgroups(int rg, std::vector<std::string> & names);
};

#endif /* TREADGROUPS_H_ */
