/*
 * TReadGroups.cpp
 *
 *  Created on: Oct 17, 2019
 *      Author: linkv
 */


#include "TReadGroups.h"

namespace BAM{


//---------------------------------------------------------------
//TReadGroup
//---------------------------------------------------------------
TReadGroup::TReadGroup(const uint16_t ID, const std::string Name){
	id = ID;
	name = Name;
	inUse = true;
	writeToHeader = true;
};

TReadGroup::TReadGroup(const TReadGroup & other){
	id = other.id;
	description = other.description;
	flowOrder = other.flowOrder;
	name = other.name;
	keySequence = other.keySequence;
	library = other.library;
	platformUnit = other.platformUnit;
	predictedInsertSize = other.predictedInsertSize;
	productionDate = other.productionDate;
	program = other.program;
	sample = other.sample;
	sequencingCenter = other.sequencingCenter;
	sequencingTechnology = other.sequencingTechnology;

	inUse = true;
	writeToHeader = true;
};

bool TReadGroup::operator<(const TReadGroup & right){
	return name < right.name;
};

bool operator<(const std::string & left, const TReadGroup & right){
	return left < right.name;
};

bool operator<(const TReadGroup & left, const std::string & right){
	return left.name < right;
};

//---------------------------------------------------------------
//TReadGroups
//---------------------------------------------------------------
TReadGroups::TReadGroups(){
	_limitReadGroups = false;
};

void TReadGroups::clear(){
	_readGroups.clear();
	_readGroupsById.clear();
};

void TReadGroups::_fillLookupFromId(){
	_readGroupsById.resize(_readGroups.size());
	for(auto& rg : _readGroups){
		_readGroupsById[rg.id] = &rg;
	}
};

TReadGroup& TReadGroups::add(const std::string name){
	auto rg = _readGroups.emplace(_readGroups.size(), name);
	if(rg.second){
		_fillLookupFromId();
	}
	return rg.first;
};

TReadGroup& TReadGroups::addAlternativeRG(const std::string Name, const std::string original){
	//getId original
	TReadGroup& rg = getReadGroup(original);

	//make sure new name does not yet exist
	if(readGroupExists(Name)){
		throw "Can not add truncated or merged read group '" + Name + "': read group already exists!";
	}

	//make copy
	TReadGroup newRg(rg);

	//set name and give new id
	newRg.name = Name;
	newRg.id = _readGroups.size();

	//add to set and inUse
	auto r = _readGroups.insert(newRg);
	if(r.second){
		_fillLookupFromId();
	}

	return r.first;
};

uint16_t TReadGroups::size() const{
	return _readGroups.size();
};

const std::string& TReadGroups::getName(uint16_t readGroupId) const{
	if(readGroupId < 0 || readGroupId >= _readGroups.size()) throw "No read group with number " + toString(readGroupId) + "!";

	return _readGroupsById[readGroupId]->name;
};

uint16_t TReadGroups::getId(const std::string & name) const{
	auto rg = _readGroups.find(name);
	if(rg != _readGroups.end())
		return rg->id;
	throw "Read Group '" + name + "' is not present in header of bam file!";
};

TReadGroup& TReadGroups::getReadGroup(const std::string & name){
	auto rg = _readGroups.find(name);
	if(rg != _readGroups.end())
		return rg;
	throw "Read Group '" + name + "' is not present in header of bam file!";
};

bool TReadGroups::readGroupExists(const std::string & name) const{
	auto rg = _readGroups.find(name);
	if(rg != _readGroups.end())
			return true;
	return false;
};

bool TReadGroups::readGroupInUse(const uint16_t & readGroupId) const{
	return _readGroupsById[readGroupId]->inUse;
};

bool TReadGroups::readGroupInUse(const std::string name) const{
	auto rg = _readGroups.find(name);
	if(rg != _readGroups.end())
		return rg->inUse;
	throw "Read Group '" + name + "' is not present in header of bam file!";
};

void TReadGroups::filterReadGroups(std::string readGroupList){
	_limitReadGroups = true;
	std::vector<std::string> readGroupsInUse;
	fillVectorFromString(readGroupList, readGroupsInUse, ',');

	//set all to false
	for(auto& rg : _readGroups){
		rg.inUse = false;
		rg.writeToHeader = false;
	}

	//set those in list to true
	for(auto& r : readGroupsInUse){
		auto& rg = _readGroups.find(r);
		if(rg == _readGroups.end())
			throw "Read Group '" + r + "' is not present in header of bam file!";
		rg->inUse = true;
		rg->writeToHeader = true;
	}
};

void TReadGroups::removeFromHeader(const std::string name){
	auto rg = _readGroups.find(name);
	if(rg == _readGroups.end())
		throw "Read Group '" + name + "' is not present in header of bam file!";
	rg->writeToHeader = false;
};

void TReadGroups::removeFromHeader(const uint16_t readGroupId){
	if(readGroupId < 0 || readGroupId >= _readGroups.size()) throw "No read group with number " + toString(readGroupId) + "!";
	_readGroupsById[readGroupId]->writeToHeader = false;
};

void TReadGroups::printReadgroupsInUse(TLog* logfile) const{
	for(auto& rg : _readGroups){
		if(rg.inUse)
			logfile->list(rg.name);
	}
};

void TReadGroups::fillVectorWithNames(std::vector<std::string> & vec) const{
	vec.resize(_readGroups.size());
	for(auto& rg : _readGroups){
		vec[rg.id] = rg.name;
	}
};

//---------------------------------------------------------------
//TReadGroupMap
//---------------------------------------------------------------
TReadGroupMap::TReadGroupMap(TReadGroups* ReadGroups){
	//no pooling: internal index = read group index
	_readGroups = ReadGroups;
	_origNumReadGroups = _readGroups->size();
	_readGroupMap = new uint16_t[_origNumReadGroups];

	_fillWithoutPooling();
};

TReadGroupMap::TReadGroupMap(TReadGroups* ReadGroups, const std::string filename, TLog* logfile){
	_readGroups = ReadGroups;
	_origNumReadGroups = _readGroups->size();
	_readGroupMap = new uint16_t[_origNumReadGroups];

	if(filename.empty()){
		_fillWithoutPooling();
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
	uint16_t oldId;

	for(uint16_t rg = 0; rg < readGroupsToMerge.size(); ++rg, ++mergeIt){
		logfile->startIndent("The following read groups will be combined into one group for parameter estimation:");
		for(std::vector<std::string>::iterator it = mergeIt->begin(); it != mergeIt->end(); ++it){
			logfile->list(*it);
			oldId = _readGroups->getId(*it);
			if(_readGroupMap[oldId] >= 0) throw "Read group '" + *it + "' is listed multiple times in file '" + filename + "'!";
			_readGroupMap[oldId] = rg;
		}
		logfile->endIndent();
	}

	_numReadGroups = readGroupsToMerge.size();

	//now add read groups that will not be merged
	bool printed = false;
	std::string name;
	for(uint16_t i = 0; i < _readGroups->size(); ++i){
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
	_reverseReadGroupMap = new std::vector<uint16_t>[_numReadGroups];
	for(uint16_t i = 0; i < _origNumReadGroups; ++i){
		_reverseReadGroupMap[getIndex(i)].push_back(i);
	}
};

uint16_t TReadGroupMap::getOrigNumReadGroups(){ return _origNumReadGroups; };
uint16_t TReadGroupMap::getNumReadGroups(){ return _numReadGroups; };

uint16_t TReadGroupMap::getIndex(const uint16_t rg){ return _readGroupMap[rg]; };
uint16_t TReadGroupMap::operator[](const uint16_t rg){ return _readGroupMap[rg]; };
uint16_t TReadGroupMap::getIndex(const std::string readGroupName){
	uint16_t rg = _readGroups->getId(readGroupName);
	return getIndex(rg);
};

void TReadGroupMap::fillNamesOfReadgroups(uint16_t rg, std::vector<std::string> & names){
	for(size_t i : _reverseReadGroupMap[rg]){
		names.push_back(_readGroups->getName(i));
	}
};

}; //end namespace
