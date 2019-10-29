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
	readGroup* groups;
	int numGroups;
	bool initialized;
	bool limitReadGroups;
	bool* inUse;

public:
	TReadGroups();
	~TReadGroups();

	void fill(BamTools::SamHeader & bamHeader);
	int find(std::string & name);

	int find(BamTools::BamAlignment & alignment);
	bool readGroupExists(std::string & name);
	bool readGroupInUse(int & readGroupId);
	bool readGroupInUse(BamTools::BamAlignment & alignment);
	std::string getName(int readGroupId);
	unsigned int size();
	void filterReadGroups(std::string readGroupList);
	void printReadgroupsInUse(TLog* logfile);
	int addTruncatedOrMergedRG(BamTools::SamHeader & bamHeader, std::string oldReadGroupName, std::string newReadGroupName);
};

//---------------------------------------------------------------
//TReadGroupMap
//---------------------------------------------------------------
class TReadGroupMap{
private:
	void initializeFromFile(TReadGroups &readGroups, std::string filename, TLog* logfile);

public:
	int origNumReadGroups;
	int numReadGroups;
	bool mergedInd;
	int* readGroupMap;

	TReadGroupMap(TReadGroups & readGroups);
	TReadGroupMap(TReadGroups & readGroups, const std::string filename, TLog* logfile);

	~TReadGroupMap();

	int getOrigNumReadGroups();
	int getNumReadGroups();

	int operator[](int rg);
	int getIndex(int rg);
};

#endif /* TREADGROUPS_H_ */
