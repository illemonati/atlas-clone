/*
 * TReadGroups.cpp
 *
 *  Created on: Oct 17, 2019
 *      Author: linkv
 */


#include "TReadGroups.h"

namespace BAM{

using coretools::str::toString;

//---------------------------------------------------------------
//TReadGroup
//---------------------------------------------------------------
TReadGroup::TReadGroup(const uint16_t ID, const std::string Name){
	id = ID;
	name_ID = Name;
	inUse = true;
	writeToHeader = true;
};

/*
TReadGroup::TReadGroup(const TReadGroup & other){
	id = other.id;
	description_DS = other.description_DS;
	flowOrder_FO = other.flowOrder_FO;
	name_ID = other.name_ID;
	keySequence_KS = other.keySequence_KS;
	library_LB = other.library_LB;
	platformUnit_PU = other.platformUnit_PU;
	predictedInsertSize_PI = other.predictedInsertSize_PI;
	productionDate_DT = other.productionDate_DT;
	program_PG = other.program_PG;
	sample_SM = other.sample_SM;
	sequencingCenter_CN = other.sequencingCenter_CN;
	sequencingTechnology_PL = other.sequencingTechnology_PL;

	inUse = true;
	writeToHeader = true;
};
*/

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

bool TReadGroup::operator<(const std::string & right) const{
	return name_ID < right;
};

bool operator<(const std::string & left, const TReadGroup & right){
	return left < right.name_ID;
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
	for(const TReadGroup& rg : _readGroups){
		int id = rg.id;
		_readGroupsById[id] = &rg;
	}
};

const TReadGroup& TReadGroups::add(const std::string name){
	const auto rg = _readGroups.emplace(_readGroups.size(), name);
	if(rg.second){
		_fillLookupFromId();
	}
	return *rg.first;
};

const TReadGroup& TReadGroups::addAlternativeRG(const std::string Name, const std::string original){
	//getId original
	const TReadGroup& rg = getReadGroup(original);

	//make sure new name does not yet exist
	if(readGroupExists(Name)){
		throw "Can not add truncated or merged read group '" + Name + "': read group already exists!";
	}

	//make copy
	TReadGroup newRg(rg);

	//set name and give new id
	newRg.name_ID = Name;
	newRg.id = _readGroups.size();

	//add to set and inUse
	auto r = _readGroups.insert(newRg);
	if(r.second){
		_fillLookupFromId();
	}

	return *r.first;
};

uint16_t TReadGroups::size() const{
	return _readGroups.size();
};

bool TReadGroups::empty() const{
	return _readGroups.empty();
};

const std::string& TReadGroups::getName(const uint16_t & readGroupId) const{
	if(readGroupId < 0 || readGroupId >= _readGroups.size()) throw "No read group with number " + toString(readGroupId) + "!";

	return _readGroupsById[readGroupId]->name_ID;
};

std::vector<std::string> TReadGroups::getNames(std::vector<uint16_t> & readGroupIds) const{
	std::vector<std::string> names;
	for(auto& r : readGroupIds){
		names.push_back(getName(r));
	}
	return names;
};


uint16_t TReadGroups::getId(const std::string & name) const{
	auto rg = _readGroups.find(name);
	if(rg != _readGroups.end())
		return rg->id;
	throw "Read Group '" + name + "' is not present in header of bam file!";
};

const TReadGroup& TReadGroups::getReadGroup(const std::string & name){
	auto rg = _readGroups.find(name);
	if(rg != _readGroups.end())
		return *rg;
	throw "Read Group '" + name + "' is not present in header of bam file!";
};

const TReadGroup& TReadGroups::operator[](const uint16_t & readGroupId) const{
	if(readGroupId < 0 || readGroupId >= _readGroups.size()) throw "No read group with number " + toString(readGroupId) + "!";
	return *_readGroupsById[readGroupId];
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
	coretools::str::fillContainerFromString(readGroupList, readGroupsInUse, ",");

	//set all to false
	for(auto& rg : _readGroups){
		rg.inUse = false;
		rg.writeToHeader = false;
	}

	//set those in list to true
	for(auto& r : readGroupsInUse){
		const auto& rg = _readGroups.find(r);
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

void TReadGroups::printReadgroupsInUse(coretools::TLog* logfile) const{
	for(auto& rg : _readGroups){
		if(rg.inUse)
			logfile->list(rg.name_ID);
	}
};

void TReadGroups::fillVectorWithNames(std::vector<std::string> & vec) const{
	vec.resize(_readGroups.size());
	for(auto& rg : _readGroups){
		vec[rg.id] = rg.name_ID;
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

TReadGroupMap::TReadGroupMap(const TReadGroups & ReadGroups, const std::string filename, coretools::TLog* logfile){
	if(filename.empty()){
		_fillWithoutPooling(ReadGroups);
	} else {
		_fillFromFile(ReadGroups, filename, logfile);
	}
};

void TReadGroupMap::_resize(const TReadGroups & ReadGroups){
	_readGroupMap.resize(ReadGroups.size(), ReadGroupMapNotInitializedIndex);
	_reverseReadGroupMap.resize(ReadGroups.size());
};

void TReadGroupMap::_markAsInUse(const uint16_t & index){
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

void TReadGroupMap::_fillFromFile(const TReadGroups & ReadGroups, const std::string & filename, coretools::TLog* logfile){
	//set all values to no-initialized
	_resize(ReadGroups);

	//read read groups and their expected lengths
	logfile->listFlush("Reading read groups to be pooled from file '" + filename + "' ...");
	coretools::TInputFile in(filename, {"readGroup", "poolWith"}, "/t", "//");

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
				throw "Read group '" + vec[1] + "' can not be pooled and pooled with at the same time!";
			}

			//check if rg can be pooled: allow self-pooling
			if(rg != pool){
				if(_readGroupMap[rg] == rg){
					throw "Read group '" + vec[0] + "' can not be pooled and pool with at the same time!";
				} else if(_readGroupMap[rg] != ReadGroupMapNotInitializedIndex){
					throw "Read group '" + vec[0] + "' is pooled multiple times in file '" + filename + "'!";
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
	logfile->done();

	//report
	if(pooledAtLeastOneRG){
		logfile->startIndent("The read groups will be pooled for parameter estimation as follows:");
		for(uint16_t r = 0; r < _readGroupMap.size(); ++r){
			if(_reverseReadGroupMap[r].size() > 0){
				logfile->startIndent(ReadGroups.getName(r) + ":");
				for(auto& s: _reverseReadGroupMap[r]){
					logfile->list(ReadGroups.getName(s));
				}
				logfile->endIndent();
			}
		}
	} else {
		logfile->warning("No read groups are pooled! Are you using the correct pooling file?");
	}
};


}; //end namespace
