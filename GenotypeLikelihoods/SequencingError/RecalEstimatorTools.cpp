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
void addUsed(std::vector<bool> &counts, size_t value) {
	if (counts.size() <= value) counts.resize(value + 1, false);
	counts[value] = true;
};

} // namespace

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
//--------------------------------------------------------------------
void TRecalDataTable::add(const BAM::TSequencedBase & base){
	++_counts;

	//add quality
	impl::addUsed(_positions, base.distFrom5Prime);
	impl::addUsed(_fragmentLengths, base.fragmentLength);
	impl::addUsed(_qualities, base.originalQuality_phredInt.get());
	impl::addUsed(_mappingQualities, base.mappingQuality.get());
};

//--------------------------------------------------------------------
// TRecalDataTables
//--------------------------------------------------------------------

void TRecalDataTables::initialize(const BAM::TReadGroupMap *ReadGroupMapObject) {
	_readGroupMap = ReadGroupMapObject;
	_tables.resize(_readGroupMap->size());
};

void TRecalDataTables::add(const std::vector<TSite> &sites) {
	for (const auto &s : sites) {
		_size += s.depth();
		if (s.depth() > 1) ++_N_g1;

		for (const auto &b : s) {
			_tables[_readGroupMap->pooledIndex(b.readGroupID)][b.mate()].add(b);
		}
	}
}

const TRecalDataTableOneReadGroup& TRecalDataTables::operator[](uint16_t readGroupId) const{
	return _tables[ _readGroupMap->pooledIndex(readGroupId) ];
};

} //end namespaceRecal GenotypeLikelihoods
