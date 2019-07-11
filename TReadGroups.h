/*
 * TReadGroups.h
 *
 *  Created on: Nov 27, 2015
 *      Author: wegmannd
 */

#ifndef TREADGROUPS_H_
#define TREADGROUPS_H_

#include "stringFunctions.h"
#include "IOTools/IOAbstractClasses/SamHeader.h"
#include "IOTools/IOAbstractClasses/BamAlignment.h"
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
	uint16_t truncatedReadGroupID;
	std::string truncatedReadGroup;

    TReadGroupMaxLength(int MaxLen, int TruncatedReadGroupID, std::string & TruncatedReadGroup);;
};

//---------------------------------------------------------------
//TReadGroups
//---------------------------------------------------------------
struct readGroup{
public:
	std::string name;
	uint16_t id;
    SamReadGroup* object;
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

    void fill(SamHeader & bamHeader);
    int find(std::string & name);
    int find(BamAlignment & alignment);
    bool readGroupExists(std::string & name);
    bool readGroupInUse(int & readGroupId);
    bool readGroupInUse(BamAlignment & alignment);
    std::string getName(int readGroupId);
    int size();
    void filterReadGroups(std::string readGroupList);
    void printReadgroupsInUse(TLog* logfile);
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
