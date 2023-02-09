/*
 * TReadGroups.cpp
 *
 *  Created on: Oct 17, 2019
 *      Author: linkv
 */


#include "TReadGroups.h"

#include <algorithm>
#include <memory>
#include <utility>

#include "coretools/Files/TFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Strings/stringFunctions.h"

namespace BAM{

using coretools::instances::logfile;

//---------------------------------------------------------------
//TReadGroup
//---------------------------------------------------------------
TReadGroup::TReadGroup(){
	_id = TReadGroups::noReadGroupId;
	name_ID = "No Read Group";
	inUse = false;
	writeToHeader = false;
};

TReadGroup::TReadGroup(const uint16_t ID, std::string_view Name){
	_id = ID;
	name_ID = Name;
	inUse = true;
	writeToHeader = true;
};

std::string TReadGroup::compileSamHeader() const{
	std::string header = "@RG\tID:" + name_ID;

	if(!barcodeSequence_BC.empty()){
		header += "\tBC:" + barcodeSequence_BC;
	}

	if(!sequencingCenter_CN.empty()){
		header += "\tCN:" + sequencingCenter_CN;
	}

	if(!description_DS.empty()){
		header += "\tDS:" + description_DS;
	}

	if(!productionDate_DT.empty()){
		header += "\tDT:" + productionDate_DT;
	}

	if(!flowOrder_FO.empty()){
		header += "\tFO:" + flowOrder_FO;
	}

	if(!keySequence_KS.empty()){
		header += "\tKS:" + keySequence_KS;
	}

	if(!library_LB.empty()){
		header += "\tLB:" + library_LB;
	}

	if(!program_PG.empty()){
		header += "\tPG:" + program_PG;
	}

	if(!predictedInsertSize_PI.empty()){
		header += "\tPI:" + predictedInsertSize_PI;
	}

	if(!sequencingTechnology_PL.empty()){
		header += "\tPL:" + sequencingTechnology_PL;
	}

	if(!platformModel_PM.empty()){
		header += "\tPM:" + platformModel_PM;
	}

	if(!platformUnit_PU.empty()){
		header += "\tPU:" + platformUnit_PU;
	}

	if(!sample_SM.empty()){
		header += "\tSM:" + sample_SM;
	}

	return header;
};

bool TReadGroup::operator<(const TReadGroup & right) const{
	return name_ID < right.name_ID;
};

bool TReadGroup::operator<(std::string_view right) const{
	return name_ID < right;
};

bool TReadGroup::operator==(std::string_view name) const{
	return name_ID == name;
};

bool operator<(std::string_view left, const TReadGroup & right){
	return left < right.name_ID;
};

//---------------------------------------------------------------
//TReadGroups
//---------------------------------------------------------------
TReadGroups::TReadGroups(){
	_limitReadGroups = false;
};

TReadGroups::TReadGroups(const TReadGroups && other){
	_readGroups = other._readGroups;
	_limitReadGroups = other._limitReadGroups;
	_fillLookupFromId();
}

TReadGroups::TReadGroups(const TReadGroups & other){
	_readGroups = other._readGroups;
	_limitReadGroups = other._limitReadGroups;
	_fillLookupFromId();
}

TReadGroups& TReadGroups::operator=(const TReadGroups & other){
	_readGroups = other._readGroups;
	_limitReadGroups = other._limitReadGroups;
	_fillLookupFromId();
	return *this;
}

TReadGroups& TReadGroups::operator=(const TReadGroups && other){
	_readGroups = std::move(other._readGroups);
	_limitReadGroups = std::move(other._limitReadGroups);
	_fillLookupFromId();
	return *this;
}

std::vector<TReadGroup>::iterator TReadGroups::_getReadGroup(std::string_view Name){
	auto rg = std::lower_bound(_readGroups.begin(),_readGroups.end(), Name);
	if(rg != _readGroups.end() && rg->name_ID == Name){
		return rg;
	}
	return _readGroups.end();
}

std::vector<TReadGroup>::const_iterator TReadGroups::_getReadGroup(std::string_view Name) const{
	auto rg = std::lower_bound(_readGroups.cbegin(),_readGroups.cend(), Name);
	if(rg != _readGroups.cend() && rg->name_ID == Name){
		return rg;
	}
	return _readGroups.cend();
}

void TReadGroups::_fillLookupFromId(){
	//sort by name
	std::sort(_readGroups.begin(), _readGroups.end());

	//fill in by ID look-up
	_readGroupsById.resize(_readGroups.size());
	for(size_t i = 0; i < _readGroups.size(); ++i){
		_readGroupsById[_readGroups[i].id()] = i;
	}
};

// add and remove read groups
void TReadGroups::clear(){
	_readGroups.clear();
	_readGroupsById.clear();
}

TReadGroup& TReadGroups::add(std::string_view Name){
	//only add if name does not yet exist
	auto rg = _getReadGroup(Name);
	if(rg == _readGroups.end()){
		_readGroups.emplace_back(_readGroups.size(), Name);
		_fillLookupFromId();
		return *_getReadGroup(Name);
	} else {
		return *rg;
	}
};

TReadGroup& TReadGroups::addAlternativeRG(std::string_view Name, std::string_view Original){
	//getId original
	const auto& rg = getReadGroup(Original);

	//make sure new name does not yet exist
	if(readGroupExists(Name)){
		UERROR("Can not add truncated or merged read group '", Name, "': read group already exists!");
	}

	//make copy
	TReadGroup newRg(rg);

	//set name and give new id
	newRg.name_ID = Name;
	newRg.setId(_readGroups.size());

	_readGroups.push_back(newRg);
	_fillLookupFromId();

	return *_getReadGroup(Name);
};

uint16_t TReadGroups::size() const{
	return _readGroups.size();
};

bool TReadGroups::empty() const{
	return _readGroups.empty();
};

// access read groups
uint16_t TReadGroups::getId(std::string_view Name) const {
	auto rg = _getReadGroup(Name);
	if(rg == _readGroups.end()){
		return noReadGroupId;
	} else {
		return rg->id();
	}
}

const TReadGroup& TReadGroups::getReadGroup(std::string_view Name) const {
	auto rg = _getReadGroup(Name);
	if(rg != _readGroups.end())
		return *rg;
	UERROR("Read Group '", Name, "' is not present in header of bam file!");
};

TReadGroup& TReadGroups::getReadGroup(std::string_view Name){
	auto rg = _getReadGroup(Name);
	if(rg != _readGroups.end())
		return *rg;
	UERROR("Read Group '", Name, "' is not present in header of bam file!");
}

const TReadGroup& TReadGroups::getReadGroup(uint16_t ReadGroupId) const {
	if(ReadGroupId >= _readGroups.size())
		UERROR("No read group with number ", ReadGroupId, "!");
	return _readGroups[ _readGroupsById[ReadGroupId] ];
}

TReadGroup& TReadGroups::getReadGroup(uint16_t ReadGroupId){
	if(ReadGroupId >= _readGroups.size())
		UERROR("No read group with number ", ReadGroupId, "!");
	return _readGroups[ _readGroupsById[ReadGroupId] ];
}


const TReadGroup& TReadGroups::operator[](uint16_t ReadGroupId) const{
	return _readGroups[_readGroupsById[ReadGroupId]];
}

bool TReadGroups::readGroupExists(std::string_view Name) const {
	return _getReadGroup(Name) != _readGroups.cend();
}

bool TReadGroups::readGroupExists(uint16_t readGroupId) const {
	return readGroupId < _readGroups.size();
}

//getters of specific entries
const std::string& TReadGroups::getName(uint16_t ReadGroupId) const{
	if(ReadGroupId >= _readGroups.size()) DEVERROR("No read group with number ", ReadGroupId, "!");
	return _readGroups[_readGroupsById[ReadGroupId]].name_ID;
}

std::vector<std::string> TReadGroups::getNames(std::vector<uint16_t> & ReadGroupIds) const{
	std::vector<std::string> names;
	for(auto& r : ReadGroupIds){
		names.push_back(getName(r));
	}
	return names;
}

//some additional tasks
void TReadGroups::filterReadGroups(std::string_view ReadGroupList){
	_limitReadGroups = true;
	std::vector<std::string> readGroupsInUse;
	coretools::str::fillContainerFromString(ReadGroupList, readGroupsInUse, ",");

	//set all to false
	for(auto& rg : _readGroups){
		rg.inUse = false;
		rg.writeToHeader = false;
	}

	//set those in list to true
	for(auto& r : readGroupsInUse){
		TReadGroup& rg = getReadGroup(r);
		rg.inUse = true;
		rg.writeToHeader = true;
	}
}

void TReadGroups::printReadgroupsInUse() const{
	for(auto& rg : _readGroups){
		if(rg.inUse)
			logfile().list(rg.name_ID);
	}
}

void TReadGroups::fillVectorWithNames(std::vector<std::string> & Vec) const{
	Vec.resize(_readGroups.size());
	for(auto& rg : _readGroups){
		Vec[rg.id()] = rg.name_ID;
	}
};

std::string TReadGroups::compileSamHeader() const{
	std::string header;
	for(auto& rg : _readGroups){
		header += rg.compileSamHeader() + "\n";
	}
	return header;
};

//---------------------------------------------------------------
//TReadGroupMap
//---------------------------------------------------------------
const uint16_t TReadGroupMap::ReadGroupMapNotInitializedIndex = 65535; //largest possible values

TReadGroupMap::TReadGroupMap(const TReadGroups & ReadGroups){
	_fillWithoutPooling(ReadGroups);
};

TReadGroupMap::TReadGroupMap(const TReadGroups & ReadGroups, std::string_view filename){
	if(filename.empty()){
		_fillWithoutPooling(ReadGroups);
	} else {
		_fillFromFile(ReadGroups, filename);
	}
};

void TReadGroupMap::_resize(const TReadGroups & ReadGroups){
	_readGroupMap.resize(ReadGroups.size(), ReadGroupMapNotInitializedIndex);
	_reverseReadGroupMap.resize(ReadGroups.size());
};

void TReadGroupMap::_markAsInUse(uint16_t index){
	_readGroupMap[index] = index;
	_reverseReadGroupMap[index].push_back(index);
	_readGroupsInUse.push_back(index);
};

void TReadGroupMap::_fillWithoutPooling(const TReadGroups & ReadGroups){
	_resize(ReadGroups)	;
	for(int r = 0; r < ReadGroups.size(); ++r){
		_markAsInUse(r);
	}
};

void TReadGroupMap::_fillFromFile(const TReadGroups & ReadGroups, std::string_view filename){
	//set all values to no-initialized
	_resize(ReadGroups);

	//read read groups and their expected lengths
	logfile().listFlush("Reading read groups to be pooled from file '", filename, "' ...");
	coretools::TInputFile in(filename, {"readGroup", "poolWith"}, "\t", "//");

	std::vector< std::vector<std::string> > readGroupsToMerge;

	//parse file and fill vectors
	std::vector<std::string> vec;
	std::string readGroup;
	bool pooledAtLeastOneRG = false;

	while(in.read(vec)){
		//ignore if read group does not exist
		if(ReadGroups.readGroupExists(vec[0])){
			//get read group index
			uint16_t rg = ReadGroups.getId(vec[0]);
			uint16_t pool = ReadGroups.getId(vec[1]);

			//check if rg to pool with is pooled itself
			if(_readGroupMap[pool] == ReadGroupMapNotInitializedIndex){
				_markAsInUse(pool);
			} else if(_readGroupMap[pool] != pool){
				UERROR("Read group '", vec[1], "' can not be pooled and pooled with at the same time!");
			}

			//check if rg can be pooled: allow self-pooling
			if(rg != pool){
				if(_readGroupMap[rg] == rg){
					UERROR("Read group '", vec[0], "' can not be pooled and pool with at the same time!");
				} else if(_readGroupMap[rg] != ReadGroupMapNotInitializedIndex){
					UERROR("Read group '", vec[0], "' is pooled multiple times in file '", filename, "'!");
				}

				//pool!
				_readGroupMap[rg] = pool;
				_reverseReadGroupMap[pool].push_back(rg);
				pooledAtLeastOneRG = true;
			}
		}
	}

	//mark all read groups not in file as pooled with itself
	for(uint16_t r = 0; r < _readGroupMap.size(); ++r){
		if(_readGroupMap[r] == ReadGroupMapNotInitializedIndex){
			_markAsInUse(r);
		}
	}
	logfile().done();

	//report
	if(pooledAtLeastOneRG){
		logfile().startIndent("The read groups will be pooled for parameter estimation as follows:");
		for(uint16_t r = 0; r < _readGroupMap.size(); ++r){
			if(_reverseReadGroupMap[r].size() > 0){
				logfile().startIndent(ReadGroups.getName(r) + ":");
				for(auto& s: _reverseReadGroupMap[r]){
					logfile().list(ReadGroups.getName(s));
				}
				logfile().endIndent();
			}
		}
	} else {
		logfile().warning("No read groups are pooled! Are you using the correct pooling file?");
	}
};


}; //end namespace
