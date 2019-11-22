/*
 * TReadGroups.cpp
 *
 *  Created on: Oct 17, 2019
 *      Author: linkv
 */


#include "TReadGroups.h"



//---------------------------------------------------------------
//TReadGroups
//---------------------------------------------------------------

TReadGroups::TReadGroups(){
	initialized = false;
	numGroups = 0;
	groups = NULL;
	inUse = NULL;
	limitReadGroups = false;
};

TReadGroups::~TReadGroups(){
	if(initialized){
		delete[] groups;
		delete[] inUse;
	}
};

void TReadGroups::fill(BamTools::SamHeader & bamHeader){
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

int TReadGroups::find(std::string & name){
	for(int i=0; i<numGroups; ++i){
		if(groups[i].name == name) return i;
	}
	throw "Read Group '" + name + "' was not present in header of bam file!";
};

int TReadGroups::find(BamTools::BamAlignment & alignment){
	std::string tmp;
	alignment.GetTag("RG", tmp);
	return find(tmp);
};

bool TReadGroups::readGroupExists(std::string & name){
	for(int i=0; i<numGroups; ++i){
		if(groups[i].name == name) return true;
	}
	return false;
};

bool TReadGroups::readGroupInUse(const int & readGroupId){
	return inUse[readGroupId];
};

bool TReadGroups::readGroupInUse(const size_t & readGroupId){
	return inUse[readGroupId];
};

bool TReadGroups::readGroupInUse(BamTools::BamAlignment & alignment){
	return inUse[find(alignment)];
};

std::string TReadGroups::getName(int readGroupId){
	if(readGroupId < 0 || readGroupId >= numGroups) throw "No read group with number " + toString(readGroupId) + "!";
	return groups[readGroupId].name;
};

unsigned int TReadGroups::size(){
	return numGroups;
};

void TReadGroups::filterReadGroups(std::string readGroupList){
	limitReadGroups = true;
	std::vector<std::string> readGroupsInUse;
	fillVectorFromString(readGroupList, readGroupsInUse, ',');
	for(int i=0; i < numGroups; i++){
		if(std::find(readGroupsInUse.begin(), readGroupsInUse.end(), getName(i)) != readGroupsInUse.end()){
			inUse[i] = true;
		} else inUse[i] = false;
	}
};

void TReadGroups::printReadgroupsInUse(TLog* logfile){
	for(int i=0; i < numGroups; i++){
		if(inUse[i])
			logfile->list(groups[i].name);
	}
};

int TReadGroups::addTruncatedOrMergedRG(BamTools::SamHeader & bamHeader, std::string origReadGroupName, std::string newReadGroupName){
	//add a new readgroup for the truncated reads to the header
	bamHeader.ReadGroups.Add(newReadGroupName);
	fill(bamHeader);
	int origReadGroupId = find(origReadGroupName);
	int newReadGroupId = find(newReadGroupName);

	//copy original tags to truncated read groups

	BamTools::SamReadGroupIterator newRG = bamHeader.ReadGroups.Begin() + newReadGroupId;
	BamTools::SamReadGroupIterator origRG = bamHeader.ReadGroups.Begin() + origReadGroupId;
	newRG->Library = origRG->Library;
	newRG->PlatformUnit = origRG->PlatformUnit;
	newRG->PredictedInsertSize = origRG->PredictedInsertSize;
	newRG->ProductionDate = origRG->ProductionDate;
	newRG->Program = origRG->Program;
	newRG->Sample = origRG->Sample;
	newRG->SequencingCenter = origRG->SequencingCenter;
	newRG->SequencingTechnology = origRG->SequencingTechnology;

	return(newReadGroupId);
};


//---------------------------------------------------------------
//TReadGroupMap
//---------------------------------------------------------------
void TReadGroupMap::initializeFromFile(TReadGroups &readGroups, std::string filename, TLog* logfile){
	//initialize to -1
	for(int i = 0; i < origNumReadGroups; ++i){
		readGroupMap[i] = -1;
	}

	//read read groups and their expected lengths
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
	logfile->done();

	std::vector< std::vector<std::string> >::iterator mergeIt = readGroupsToMerge.begin();
	int oldId;

	for(unsigned int rg = 0; rg < readGroupsToMerge.size(); ++rg, ++mergeIt){
		logfile->startIndent("The following read groups will be combined into one group for parameter estimation:");
		for(std::vector<std::string>::iterator it = mergeIt->begin(); it != mergeIt->end(); ++it){
			logfile->list(*it);
			oldId = readGroups.find(*it);
			if(readGroupMap[oldId] >= 0) throw "Read group '" + *it + "' is listed multiple times in file '" + filename + "'!";
			readGroupMap[oldId] = rg;
		}
		logfile->endIndent();
	}

	numReadGroups = readGroupsToMerge.size();

	//now add read groups that will not be merged
	bool printed = false;
	std::string name;
	for(size_t i = 0; i < readGroups.size(); ++i){
		//check if it is mapped, otherwise add
		if(readGroupMap[i] < 0){
			if(!printed){
				logfile->startIndent("The following read groups will be kept as is:");
				printed = true;
			}
			name = readGroups.getName(i);
			logfile->list(name);
			readGroupMap[i] = numReadGroups;
			++numReadGroups;
		}
	}

	if(printed) logfile->endIndent();
	else logfile->list("All existing read groups will be merged into a new read group.");
};



TReadGroupMap::TReadGroupMap(TReadGroups & readGroups){
	origNumReadGroups = readGroups.size();
	readGroupMap = new int[origNumReadGroups];
	mergedInd = false;
	numReadGroups = origNumReadGroups;
	for(int i = 0; i < numReadGroups; ++i){
		readGroupMap[i] = i;
	}
};

TReadGroupMap::TReadGroupMap(TReadGroups & readGroups, const std::string filename, TLog* logfile){
	origNumReadGroups = readGroups.size();
	readGroupMap = new int[origNumReadGroups];
	mergedInd = false;

	if(filename.empty()){
		numReadGroups = origNumReadGroups;
		for(int i = 0; i < numReadGroups; ++i){
			readGroupMap[i] = i;
		}
	} else {
		initializeFromFile(readGroups, filename, logfile);
	}
};

TReadGroupMap::~TReadGroupMap(){
	delete[] readGroupMap;
};

int TReadGroupMap::getOrigNumReadGroups(){ return origNumReadGroups; };
int TReadGroupMap::getNumReadGroups(){ return numReadGroups; };

int TReadGroupMap::operator[](int rg){ return readGroupMap[rg]; };
int TReadGroupMap::getIndex(int rg){ return readGroupMap[rg]; };

