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
	uint16_t truncatedReadGroupID;
	std::string truncatedReadGroup;

	TReadGroupMaxLength(int MaxLen, int TruncatedReadGroupID, std::string & TruncatedReadGroup){
		maxLen = MaxLen;
		truncatedReadGroupID = TruncatedReadGroupID;
		truncatedReadGroup = TruncatedReadGroup;
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
	TReadGroups(){
		initialized = false;
		numGroups = 0;
		groups = NULL;
		inUse = NULL;
		limitReadGroups = false;
	};

	~TReadGroups(){
		if(initialized){
			delete[] groups;
			delete[] inUse;
		}
	};

	void fill(BamTools::SamHeader & bamHeader){
		//empty if filled before
		if(initialized) delete[] groups;
		//create and fill array
		numGroups = bamHeader.ReadGroups.Size();
		groups = new readGroup[numGroups];
		inUse = new bool[numGroups];
		int i = 0;
		for(BamTools::SamReadGroupIterator it = bamHeader.ReadGroups.Begin(); it != bamHeader.ReadGroups.End(); ++it, ++i){
			groups[i].id = i;
			groups[i].name = it->ID;
			groups[i].object= &(*it);
			inUse[i] = true;
		}
		initialized = true;
	};

	int find(std::string & name){
		for(int i=0; i<numGroups; ++i){
			if(groups[i].name == name) return i;
		}
		throw "Read Group '" + name + "' was not present in header of bam file!";
	};

	int find(BamTools::BamAlignment & alignment){
		std::string tmp;
		alignment.GetTag("RG", tmp);
		return find(tmp);
	};

	bool readGroupExists(std::string & name){
		for(int i=0; i<numGroups; ++i){
			if(groups[i].name == name) return true;
		}
		return false;
	};

	bool readGroupInUse(int & readGroupId){
		return inUse[readGroupId];
	};

	bool readGroupInUse(BamTools::BamAlignment & alignment){
		return inUse[find(alignment)];
	};

	std::string getName(int readGroupId){
		if(readGroupId < 0 || readGroupId >= numGroups) throw "No read group with number " + toString(readGroupId) + "!";
		return groups[readGroupId].name;
	};

	int size(){
		return numGroups;
	};

	void filterReadGroups(std::string readGroupList){
		limitReadGroups = true;
		std::vector<std::string> readGroupsInUse;
		fillVectorFromString(readGroupList, readGroupsInUse, ',');
		for(int i=0; i < numGroups; i++){
			if(std::find(readGroupsInUse.begin(), readGroupsInUse.end(), getName(i)) != readGroupsInUse.end()){
				inUse[i] = true;
			} else inUse[i] = false;
		}
	};

	void printReadgroupsInUse(TLog* logfile){
		for(int i=0; i < numGroups; i++){
			if(inUse[i])
				logfile->list(groups[i].name);
		}
	};
};


//---------------------------------------------------------------
//TReadGroupMap
//---------------------------------------------------------------
class TReadGroupMap{
private:

public:
	int origNumReadGroups;
	int numReadGroups;
	bool readGroupMapInitialized;
	bool mergedInd;
	int* readGroupMap;

	TReadGroupMap(BamTools::SamHeader* bamHeader, TParameters & params, TLog* logfile){
		origNumReadGroups = bamHeader->ReadGroups.Size();
		readGroupMap = new int[origNumReadGroups];
		readGroupMapInitialized = true;
		if(params.parameterExists("poolReadGroups")) mergedInd = true;
		else mergedInd = false;

		//construct array from vectors and report
		for(int i=0; i<origNumReadGroups; ++i)	readGroupMap[i] = -1; //map initialized

		if(mergedInd){
			//read read groups and their expected lengths
			std::string filename = params.getParameterString("poolReadGroups");
			if(filename=="") throw "No file specifying read groups to merge provided!";
			logfile->listFlush("Reading read groups to be merged from file '" + filename + "' ...");
			std::vector< std::vector<std::string> > readGroupsToMerge;
			std::vector< std::vector<std::string> >::reverse_iterator rIt;
			std::ifstream file(filename.c_str());
			if(!file) throw "Failed to open file '" + filename + "!";

			//parse file and fill vectors
			int lineNum = 0;
			std::vector<std::string> vec;
			std::string readGroup;
			while(file.good() && !file.eof()){
				++lineNum;
				fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
				if(!vec.empty()){
					if(vec.size() < 2) throw "Wrong number of entries on line " + toString(lineNum) + " in file '" + filename + "'! Read groups cannot be merged with themselves!";
					//add to new header
					//others are those to be merged: find read group in header and store int
					readGroupsToMerge.push_back(std::vector<std::string>());
					rIt = readGroupsToMerge.rbegin();
					for(unsigned int i=0; i<vec.size(); ++i){
						rIt->push_back(vec[i]);
					}
				}
			}
			TReadGroups ReadGroupObject;
			ReadGroupObject.fill(*bamHeader);
			logfile->done();

			std::vector< std::vector<std::string> >::iterator mergeIt = readGroupsToMerge.begin();
			int oldId;

			for(unsigned int rg = 0; rg < readGroupsToMerge.size(); ++rg, ++mergeIt){
				logfile->startIndent("The following read groups will be combined into one group for parameter estimation:");
				for(std::vector<std::string>::iterator it = mergeIt->begin(); it != mergeIt->end(); ++it){
					logfile->list(*it);
					oldId = ReadGroupObject.find(*it);
					if(readGroupMap[oldId] >= 0) throw "Read group '" + *it + "' is listed multiple times in file '" + filename + "'!";
					readGroupMap[oldId] = rg;
				}
				logfile->endIndent();
			}

			numReadGroups = readGroupsToMerge.size();

			//now add read groups that will not be merged
			bool printed = false;
			std::string name;
			for(int i = 0; i < ReadGroupObject.size(); ++i){
				//check if it is mapped, otherwise add
				if(readGroupMap[i] < 0){
					if(!printed){
						logfile->startIndent("The following read groups will be kept as is:");
						printed = true;
					}
					name = ReadGroupObject.getName(i);
					logfile->list(name);
					readGroupMap[i] = numReadGroups;
					++numReadGroups;
				}
			}

			if(printed) logfile->endIndent();
			else logfile->list("All existing read groups will be merged into a new read group.");
		} else {
			numReadGroups = origNumReadGroups;
			for(int i = 0; i < numReadGroups; ++i){
				readGroupMap[i] = i;
			}
		}
	}

	~TReadGroupMap(){
		if(readGroupMapInitialized) delete[] readGroupMap;
	}

	int getOrigNumReadGroups(){ return origNumReadGroups; }
	int getNumReadGroups(){ return numReadGroups; }

	int operator[](int rg){
		return readGroupMap[rg];
	};
};

#endif /* TREADGROUPS_H_ */
