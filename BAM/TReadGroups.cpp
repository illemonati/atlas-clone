/*
 * TReadGroups.cpp
 *
 *  Created on: Oct 17, 2019
 *      Author: linkv
 */


#include "TReadGroups.h"

#include "coretools/Main/TLog.h"
#include "coretools/Strings/fillContainer.h"

namespace BAM{

using coretools::instances::logfile;

//---------------------------------------------------------------
//TReadGroups
//---------------------------------------------------------------

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
		_readGroupsById[_readGroups[i].id] = i;
	}
}

TReadGroups::TReadGroups() {
	_readGroupIdForReadsWithoutReadGroup = noReadGroupId;
}

// add and remove read groups
void TReadGroups::clear(){
	_readGroups.clear();
	_readGroupsById.clear();
	_readGroupIdForReadsWithoutReadGroup = noReadGroupId;
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
}

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
	newRg.id      = _readGroups.size();

	_readGroups.push_back(newRg);
	_fillLookupFromId();

	return *_getReadGroup(Name);
}
void TReadGroups::addReadGroupForReadsWithoutReadGroup(std::string_view Name){
	if(_readGroupIdForReadsWithoutReadGroup == noReadGroupId){
		auto& rg = add(Name);
		_readGroupIdForReadsWithoutReadGroup = rg.id;
	}
}

size_t TReadGroups::size() const{
	return _readGroups.size();
}

bool TReadGroups::empty() const{
	return _readGroups.empty();
}

// access read groups
size_t TReadGroups::getId(std::string_view Name) const {
	auto rg = _getReadGroup(Name);
	if(rg == _readGroups.end()){		
		return _readGroupIdForReadsWithoutReadGroup;
	} else {
		return rg->id;
	}
}

const TReadGroup& TReadGroups::getReadGroup(std::string_view Name) const {
	auto rg = _getReadGroup(Name);
	if(rg != _readGroups.end())
		return *rg;
	UERROR("Read Group '", Name, "' is not present in header of bam file!");
}

TReadGroup& TReadGroups::getReadGroup(std::string_view Name){
	auto rg = _getReadGroup(Name);
	if(rg != _readGroups.end())
		return *rg;
	UERROR("Read Group '", Name, "' is not present in header of bam file!");
}

const TReadGroup& TReadGroups::getReadGroup(size_t ReadGroupId) const {
	if (ReadGroupId == noReadGroupId) return _noReadGroup;
	if(ReadGroupId >= _readGroups.size())
		UERROR("No read group with number ", ReadGroupId, "!");
	return _readGroups[ _readGroupsById[ReadGroupId] ];
}

TReadGroup& TReadGroups::getReadGroup(size_t ReadGroupId){
	if(ReadGroupId >= _readGroups.size())
		UERROR("No read group with number ", ReadGroupId, "!");
	return _readGroups[ _readGroupsById[ReadGroupId] ];
}


const TReadGroup& TReadGroups::operator[](size_t ReadGroupId) const{
	return _readGroups[_readGroupsById[ReadGroupId]];
}

bool TReadGroups::readGroupExists(std::string_view Name) const {
	return _getReadGroup(Name) != _readGroups.cend();
}

//getters of specific entries
const std::string& TReadGroups::getName(size_t ReadGroupId) const{
	if(ReadGroupId >= _readGroups.size()) DEVERROR("No read group with number ", ReadGroupId, "!");
	return _readGroups[_readGroupsById[ReadGroupId]].name_ID;
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
		Vec[rg.id] = rg.name_ID;
	}
}

std::string TReadGroups::compileSamHeader() const{
	std::string header;
	for(auto& rg : _readGroups){
		if (rg.writeToHeader) header += rg.compileSamHeader() + "\n";
	}
	return header;
}


}; //end namespace
