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
		t.front().clear();
		t.back().clear();
	}
	_totalCounts = 0;
};

void TRecalDataTables::add(const BAM::TSequencedBase & base){
	++_totalCounts;
	_tables[_readGroupMap->pooledIndex(base.readGroupID)][base.mate()].add(base);
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

} //end namespaceRecal GenotypeLikelihoods
