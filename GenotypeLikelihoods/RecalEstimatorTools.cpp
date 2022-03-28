/*
 * RecalEstimatorTools.cpp
 *
 *  Created on: May 4, 2021
 *      Author: phaentu
 */

#include "RecalEstimatorTools.h"
#include "TLog.h"

namespace GenotypeLikelihoods {

namespace RecalEstimatorTools {

using coretools::instances::logfile;

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
//--------------------------------------------------------------------
TRecalDataTable::TRecalDataTable(){
	_counts = 0;
};

void TRecalDataTable::add(const BAM::TSequencedBase & base){
	++_counts;

	//add quality
	_positions.add(base.distFrom5Prime);
	_fragmentLengths.add(base.fragmentLength);
	_qualities.add(base.originalQuality_phredInt.get());
	_mappingQualities.add(base.mappingQuality);
};

void TRecalDataTable::clear(){
	_counts = 0;
	_positions.clear();
	_fragmentLengths.clear();
	_qualities.clear();
	_mappingQualities.clear();
};

uint64_t TRecalDataTable::size() const{
	return _counts;
};

const TRecalDataVector<uint16_t>& TRecalDataTable::positions() const{
	return _positions;
};

const TRecalDataVector<uint16_t>& TRecalDataTable::fragmentLengths() const{
	return _fragmentLengths;
};

const TRecalDataVector<uint16_t>& TRecalDataTable::qualities() const{
	return _qualities;
};

const TRecalDataVector<uint16_t>& TRecalDataTable::mappingQualities() const{
	return _mappingQualities;
};

//--------------------------------------------------------------------
// TRecalDataTableOneReadGroup
//--------------------------------------------------------------------
const TRecalDataTable& TRecalDataTableOneReadGroup::operator[](const bool & IsSecondMate) const{
	return _tables[(int) IsSecondMate];
};

void TRecalDataTableOneReadGroup::add(const BAM::TSequencedBase & base){
	_tables[base.isSecondMate()].add(base);
};

void TRecalDataTableOneReadGroup::clear(){;
	_tables[0].clear();
	_tables[1].clear();
};

//--------------------------------------------------------------------
// TRecalDataTables
//--------------------------------------------------------------------
TRecalDataTables::TRecalDataTables(){
	_readGroups = nullptr;
	_readGroupMap = nullptr;
	_totalCounts = 0;
};

TRecalDataTables::TRecalDataTables(const BAM::TReadGroups* ReadGroups, const BAM::TReadGroupMap* ReadGroupMapObject){
	initialize(ReadGroups, ReadGroupMapObject);
};

void TRecalDataTables::initialize(const BAM::TReadGroups* ReadGroups, const BAM::TReadGroupMap* ReadGroupMapObject){
	clear();
	_readGroups = ReadGroups;
	_readGroupMap = ReadGroupMapObject;
	_tables.resize(_readGroupMap->numReadGroupsInUse());
};

void TRecalDataTables::clear(){
	for(auto& t : _tables){
		t.clear();
	}
	_totalCounts = 0;
};

void TRecalDataTables::add(const BAM::TSequencedBase & base){
	++_totalCounts;
	_tables[ _readGroupMap->pooledIndex( base.readGroupID) ].add(base);
};

void TRecalDataTables::add(const TSite & Site){
	_totalCounts += Site.depth();
	for(std::vector<BAM::TSequencedBase>::const_iterator it = Site.cbegin(); it != Site.cend(); ++it){
		add(*it);
	}
};

void TRecalDataTables::add(const std::vector<TSite> & sites){
	for(auto& s : sites){
		add(s);
	}
};

uint64_t TRecalDataTables::size() const{
	return _totalCounts;
};

const TRecalDataTableOneReadGroup& TRecalDataTables::operator[](uint16_t readGroupId) const{
	return _tables[ _readGroupMap->pooledIndex(readGroupId) ];
};

//------------------------------------------------
// Classes to keep track of models to estimate
//------------------------------------------------

std::string TModelStatusEntry::getString() const{
	if(_first && _second){
		return "(first and second mates)";
	} else if(_first){
		return "(first mate)";
	} else if(_second){
		return "(second mate)";
	} else {
		return "none";
	}
};

TModelStatusEntry& TModelStatus::operator[](const ModelStatusTypes & Type){
	return _status[static_cast<size_t>(Type)];
};

void TModelStati::add(uint16_t ReadGroupId){
	modelStatus.emplace(std::pair<uint16_t, TModelStatus>(ReadGroupId, TModelStatus()));
};

TModelStatus& TModelStati::operator[](uint16_t ReadGroupId){
	return modelStatus[ReadGroupId];
};

uint16_t TModelStati::num(const ModelStatusTypes & Type){
	uint16_t num = 0;
	for(auto& m : modelStatus){
		m.second[Type].size();
	}
	return num;
};

void TModelStati::report(const ModelStatusTypes & Type, const std::string & Title, const BAM::TReadGroups & ReadGroups){
	if(num(Type) > 0){
		logfile().startIndent(Title);
		for(auto& m : modelStatus){
			if(m.second[Type].size() > 0){
				logfile().list(ReadGroups.getName(m.first), " ", m.second[Type].getString());
			}
		}
		logfile().endIndent();
	}
};

}; //end namespaceRecal GenotypeLikelihoods

}; //end namespaceRecal EstimatorTools



