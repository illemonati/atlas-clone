/*
 * RecalEstimatorTools.cpp
 *
 *  Created on: May 4, 2021
 *      Author: phaentu
 */

#include "RecalEstimatorTools.h"
#include "TLog.h"
#include <numeric>

namespace GenotypeLikelihoods {

namespace RecalEstimatorTools {

using coretools::instances::logfile;


std::vector<uint16_t> vectorOfUsed(const std::vector<uint32_t> & counts) {
	//insert all with counts
	std::vector<uint16_t> vec;
	for(size_t i = 0; i < counts.size(); ++i){
		if(counts[i] > 0){
			vec.push_back(i);
		}
	}
	return vec;
};

namespace /*anonymous*/ {
void addCount(std::vector<uint32_t> &counts, uint16_t value) {
	if (counts.size() <= value) counts.resize(value + 1, 0);
	++counts[value];
};

} // namespace

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
//--------------------------------------------------------------------
void TRecalDataTable::add(const BAM::TSequencedBase & base){
	++_counts;

	//add quality
	addCount(_positions, base.distFrom5Prime);
	addCount(_fragmentLengths, base.fragmentLength);
	addCount(_qualities, base.originalQuality_phredInt.get());
	addCount(_mappingQualities, base.mappingQuality);
};

void TRecalDataTable::clear() noexcept {
	_counts = 0;
	_positions.clear();
	_fragmentLengths.clear();
	_qualities.clear();
	_mappingQualities.clear();
};

const std::vector<uint32_t>& TRecalDataTable::positions() const{
	return _positions;
};

const std::vector<uint32_t>& TRecalDataTable::fragmentLengths() const{
	return _fragmentLengths;
};

const std::vector<uint32_t>& TRecalDataTable::qualities() const{
	return _qualities;
};

const std::vector<uint32_t>& TRecalDataTable::mappingQualities() const{
	return _mappingQualities;
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
		t[0].clear();
		t[1].clear();
	}
	_totalCounts = 0;
};

void TRecalDataTables::add(const BAM::TSequencedBase & base){
	++_totalCounts;
	_tables[_readGroupMap->pooledIndex(base.readGroupID)][base.isSecondMate()].add(base);
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

std::string TModelStatusEntry::getString() const {
	if (!_bs.to_ulong()) return "none";

	if (_bs.get<0>()) {
		if (_bs.get<1>())
			return "(first and second mates)";
		else
			return "(first mate)";
	} else {
		return "(second mate)";
	}
};

void TModelStati::add(uint16_t ReadGroupId){
	modelStatus.emplace(std::pair<uint16_t, TModelStatus>(ReadGroupId, TModelStatus()));
};

TModelStatus& TModelStati::operator[](uint16_t ReadGroupId){
	return modelStatus[ReadGroupId];
};

uint16_t TModelStati::num(ModelStatusTypes Type) {
	return std::accumulate(modelStatus.cbegin(), modelStatus.cend(), 0,
			       [Type](auto tot, auto p) { return tot + p.second[Type].size(); });
};

void TModelStati::report(ModelStatusTypes Type, const std::string & Title, const BAM::TReadGroups & ReadGroups){
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



