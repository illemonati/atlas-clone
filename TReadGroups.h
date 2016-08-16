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

//---------------------------------------------------------------
//TReadGroupMaxLength
//---------------------------------------------------------------
struct TReadGroupMaxLength{
public:
	int maxLen;
	int truncatedReadGroupID;
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
	int id;
	BamTools::SamReadGroup* object;
};

//---------------------------------------------------------------
//TReadGroups
//---------------------------------------------------------------

class TReadGroups{
public:
	readGroup* groups;
	int numGroups;
	bool initialized;

	TReadGroups(){
		initialized = false;
		numGroups = 0;
		groups = NULL;
	};

	~TReadGroups(){
		if(initialized) delete[] groups;
	};

	void fill(BamTools::SamHeader & bamHeader){
		//empty if filled before
		if(initialized) delete[] groups;
		//create and fill array
		numGroups = bamHeader.ReadGroups.Size();
		groups = new readGroup[numGroups];
		int i = 0;
		for(BamTools::SamReadGroupIterator it = bamHeader.ReadGroups.Begin(); it != bamHeader.ReadGroups.End(); ++it, ++i){
			groups[i].id = i;
			groups[i].name = it->ID;
			groups[i].object= &(*it);
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

	std::string getName(int num){
		if(num < 0 || num >= numGroups) throw "No read group with number " + toString(num) + "!";
		return groups[num].name;
	};

	int size(){
		return numGroups;
	};
};

#endif /* TREADGROUPS_H_ */
