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
	_initialized = false;
	_numGroups = 0;
	_groups = NULL;
	_inUse = NULL;
	_limitReadGroups = false;
};

TReadGroups::~TReadGroups(){
	if(_initialized){
		delete[] _groups;
		delete[] _inUse;
	}
};

void TReadGroups::fill(BamTools::SamHeader & bamHeader){
	//empty if filled before
	if(_initialized) delete[] _groups;
	//create and fill array
	_numGroups = bamHeader.ReadGroups.Size();
	_groups = new readGroup[_numGroups];
	_inUse = new bool[_numGroups];
	int i = 0;
	for(BamTools::SamReadGroupIterator it = bamHeader.ReadGroups.Begin(); it != bamHeader.ReadGroups.End(); ++it, ++i){
		_groups[i].id = i;
		_groups[i].name = it->ID;
		_groups[i].object= &(*it);
		_inUse[i] = true;
	}
	_initialized = true;
};

int TReadGroups::find(const std::string & name){
	for(size_t i=0; i<_numGroups; ++i){
		if(_groups[i].name == name) return i;
	}
	throw "Read Group '" + name + "' was not present in header of bam file!";
};

int TReadGroups::find(BamTools::BamAlignment & alignment){
	std::string tmp;
	alignment.GetTag("RG", tmp);
	return find(tmp);
};

bool TReadGroups::readGroupExists(const std::string & name){
	for(size_t i=0; i<_numGroups; ++i){
		if(_groups[i].name == name) return true;
	}
	return false;
};

bool TReadGroups::readGroupInUse(const int & readGroupId){
	return _inUse[readGroupId];
};

bool TReadGroups::readGroupInUse(const size_t & readGroupId){
	return _inUse[readGroupId];
};

bool TReadGroups::readGroupInUse(const std::string name){
	if(!readGroupExists(name))
		return false;
	int readGroupId = find(name);
	return _inUse[readGroupId];
};

bool TReadGroups::readGroupInUse(BamTools::BamAlignment & alignment){
	return _inUse[find(alignment)];
};

std::string TReadGroups::getName(int readGroupId){
	if(readGroupId < 0 || (size_t) readGroupId >= _numGroups) throw "No read group with number " + toString(readGroupId) + "!";
	return _groups[readGroupId].name;
};

size_t TReadGroups::size(){
	return _numGroups;
};

void TReadGroups::filterReadGroups(std::string readGroupList){
	_limitReadGroups = true;
	std::vector<std::string> readGroupsInUse;
	fillVectorFromString(readGroupList, readGroupsInUse, ',');
	for(size_t i=0; i < _numGroups; i++){
		if(std::find(readGroupsInUse.begin(), readGroupsInUse.end(), getName(i)) != readGroupsInUse.end()){
			_inUse[i] = true;
		} else _inUse[i] = false;
	}
};

void TReadGroups::printReadgroupsInUse(TLog* logfile){
	for(size_t i=0; i < _numGroups; i++){
		if(_inUse[i])
			logfile->list(_groups[i].name);
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
TReadGroupMap::TReadGroupMap(TReadGroups* ReadGroups){
	//no pooling: internal index = read group index
	_readGroups = ReadGroups;
	_origNumReadGroups = _readGroups->size();

	_fillWithoutPooling();
};

TReadGroupMap::TReadGroupMap(TReadGroups* ReadGroups, const std::string filename, TLog* logfile){
	_readGroups = ReadGroups;
	_origNumReadGroups = _readGroups->size();
	_readGroupMap = new int[_origNumReadGroups];

	if(filename.empty()){
		//no pooling: internal index = read group index
		_numReadGroups = _origNumReadGroups;
		for(int i = 0; i < _numReadGroups; ++i){
			_readGroupMap[i] = i;
		}
	} else {
		_fillFromFile(filename, logfile);
	}
};

TReadGroupMap::~TReadGroupMap(){
	delete[] _readGroupMap;
	delete[] _reverseReadGroupMap;
};

void TReadGroupMap::_fillWithoutPooling(){
	_numReadGroups = _origNumReadGroups;
	_readGroupMap = new int[_origNumReadGroups];
	for(int i = 0; i < _numReadGroups; ++i){
		_readGroupMap[i] = i;
	}

	//fill reverse map
	_fillReverseMap();
};

void TReadGroupMap::_fillFromFile(std::string filename, TLog* logfile){
	//initialize to -1
	for(int i = 0; i < _origNumReadGroups; ++i){
		_readGroupMap[i] = -1;
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
			oldId = _readGroups->find(*it);
			if(_readGroupMap[oldId] >= 0) throw "Read group '" + *it + "' is listed multiple times in file '" + filename + "'!";
			_readGroupMap[oldId] = rg;
		}
		logfile->endIndent();
	}

	_numReadGroups = readGroupsToMerge.size();

	//now add read groups that will not be merged
	bool printed = false;
	std::string name;
	for(size_t i = 0; i < _readGroups->size(); ++i){
		//check if it is mapped, otherwise add
		if(_readGroupMap[i] < 0){
			if(!printed){
				logfile->startIndent("The following read groups will be kept as is:");
				printed = true;
			}
			name = _readGroups->getName(i);
			logfile->list(name);
			_readGroupMap[i] = _numReadGroups;
			++_numReadGroups;
		}
	}

	if(printed) logfile->endIndent();
	else logfile->list("All existing read groups will be merged into a new read group.");

	//fill reverse map
	_fillReverseMap();
};

void TReadGroupMap::_fillReverseMap(){
	_reverseReadGroupMap = new std::vector<int>[_numReadGroups];
	for(int i = 0; i < _origNumReadGroups; ++i){
		_reverseReadGroupMap[getIndex(i)].push_back(i);
	}
};

int TReadGroupMap::getOrigNumReadGroups(){ return _origNumReadGroups; };
int TReadGroupMap::getNumReadGroups(){ return _numReadGroups; };

int TReadGroupMap::getIndex(int rg){ return _readGroupMap[rg]; };
int TReadGroupMap::getIndex(const std::string readGroupName){
	int rg = _readGroups->find(readGroupName);
	return getIndex(rg);
};

void TReadGroupMap::fillNamesOfReadgroups(int rg, std::vector<std::string> & names){
	for(int i : _reverseReadGroupMap[rg]){
		names.push_back(_readGroups->getName(i));
	}
};

