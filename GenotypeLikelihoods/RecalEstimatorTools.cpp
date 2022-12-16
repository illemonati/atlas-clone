/*
 * RecalEstimatorTools.cpp
 *
 *  Created on: May 4, 2021
 *      Author: phaentu
 */

#include "RecalEstimatorTools.h"

#include <algorithm>
#include <numeric>
#include <utility>

#include "genometools/PhredProbabilityTypes.h"
#include "coretools/Main/TLog.h"
#include "TSequencedBase.h"
#include "TSite.h"

namespace GenotypeLikelihoods::RecalEstimatorTools {

using coretools::instances::logfile;

std::vector<uint16_t> indexedRange(const std::vector<uint32_t> & counts) {
	//insert all with counts
	std::vector<uint16_t> vec;
	for(size_t i = 0; i < counts.size(); ++i){
		if(counts[i] > 0){
			vec.push_back(i);
		}
	}
	return vec;
};

std::vector<uint16_t> fullRange(const std::vector<uint32_t>& counts) {
		const auto N = counts.size();
		std::vector<uint16_t> v(N);
		std::iota(v.begin(), v.end(), uint16_t{});
		return v;
}

namespace impl {
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
	impl::addCount(_positions, base.distFrom5Prime);
	impl::addCount(_fragmentLengths, base.fragmentLength);
	impl::addCount(_qualities, base.originalQuality_phredInt.get());
	impl::addCount(_mappingQualities, base.mappingQuality.get());
};

void TRecalDataTable::clear() noexcept {
	_counts = 0;
	_positions.clear();
	_fragmentLengths.clear();
	_qualities.clear();
	_mappingQualities.clear();
};

//--------------------------------------------------------------------
// TRecalDataTables
//--------------------------------------------------------------------

void TRecalDataTables::initialize(const BAM::TReadGroups *ReadGroups, const BAM::TReadGroupMap *ReadGroupMapObject) {
	clear();
	_readGroups   = ReadGroups;
	_readGroupMap = ReadGroupMapObject;
	_tables.resize(_readGroups->size());
};

void TRecalDataTables::clear() {
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

void TRecalDataTables::add(const TSite &Site) {
	for (const auto &b : Site) add(b);
};

void TRecalDataTables::add(const std::vector<TSite> &sites) {
	for (const auto &s : sites) add(s);
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
	modelStatus.emplace(ReadGroupId, TModelStatus());
};

TModelStatus& TModelStati::operator[](uint16_t ReadGroupId){
	return modelStatus[ReadGroupId];
};

uint16_t TModelStati::num(ModelStatusTypes Type) const {
	return std::accumulate(modelStatus.cbegin(), modelStatus.cend(), 0,
			       [Type](auto tot, auto p) { return tot + p.second[Type].size(); });
};

void TModelStati::report(ModelStatusTypes Type, const std::string & Title, const BAM::TReadGroups & ReadGroups) const {
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

} //end namespaceRecal GenotypeLikelihoods
