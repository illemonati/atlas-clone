/*
 * RecalEstimatorTools.cpp
 *
 *  Created on: May 4, 2021
 *      Author: phaentu
 */

#include "RecalEstimatorTools.h"
#include "TSite.h"
#include "coretools/Files/TOutputFile.h"

namespace GenotypeLikelihoods::RecalEstimatorTools {

namespace impl {
void addUsed(std::vector<size_t> &counts, size_t value) {
	if (counts.size() <= value) counts.resize(value + 1, 0);
	++counts[value];
}

} // namespace impl

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
//--------------------------------------------------------------------
void TRecalDataTable::add(const BAM::TSequencedBase & base){
	++_counts;

	//add quality
	impl::addUsed(_tables[SequencingError::Covariates::Context], coretools::index(base.previousBase));
	impl::addUsed(_tables[SequencingError::Covariates::FragmentLength], base.fragmentLength);
	impl::addUsed(_tables[SequencingError::Covariates::MappingQuality], base.mappingQuality.get());
	impl::addUsed(_tables[SequencingError::Covariates::Position], base.distFrom5Prime);
	impl::addUsed(_tables[SequencingError::Covariates::Quality], base.originalQuality.get());
}

//--------------------------------------------------------------------
// TRecalDataTables
//--------------------------------------------------------------------

void TRecalDataTables::add(const TSite &site) {
	_size += site.depth();
	if (site.depth() > 1) ++_N_g1;

	for (const auto &b : site) { _tables[_readGroupMap->pooledIndex(b.readGroupID)][b.mate()].add(b); }
}

const TRecalDataTableOneReadGroup& TRecalDataTables::operator[](size_t readGroupId) const{
	return _tables[ _readGroupMap->pooledIndex(readGroupId) ];
}

void TRecalDataTables::write(std::string_view Name) const {
	using coretools::TStrongArray;
	using SequencingError::Covariates;
	TStrongArray<coretools::TOutputFile, Covariates> files;

	files[Covariates::Context].open(std::string(Name).append("_contexts.txt.gz"));
	files[Covariates::Context].write("context");

	files[Covariates::FragmentLength].open(std::string(Name).append("_fragmentLengths.txt.gz"));
	files[Covariates::FragmentLength].write("fragmentLength");

	files[Covariates::MappingQuality].open(std::string(Name).append("_mappingQualities.txt.gz"));
	files[Covariates::MappingQuality].write("mappingQuality");

	files[Covariates::Position].open(std::string(Name).append("_positions.txt.gz"));
	files[Covariates::Position].write("position");

	files[Covariates::Quality].open(std::string(Name).append("_qualities.txt.gz"));
	files[Covariates::Quality].write("quality");

	std::vector<const TRecalDataTable*> usedTables;
	for (auto rg: _readGroupMap->readGroupsInUse()) {
		const auto & table = _tables[rg];
		if (table.front().size() > 0) {
			for (auto& f: files) f.writeNoDelim("RG_", rg, "_0").writeDelim();
			usedTables.push_back(&table.front());
		}
		if (table.back().size() > 0) {
			for (auto& f: files) f.writeNoDelim("RG_", rg, "_0").writeDelim();
			usedTables.push_back(&table.back());
		}
	}
	for (auto& f: files) f.endln();

	for (auto c = Covariates::min; c < Covariates::max; ++c) {
		size_t N = 0;
		for (auto pt : usedTables) {
			const auto & table = *pt;
			N = std::max(N, table[c].size());
		}
		for (size_t i = 0; i < N; ++i) {
			files[c].write(i);
			for (auto pt : usedTables) {
				const auto & table = *pt;
				if (i < table[c].size()) files[c].write(table[c][i]);
				else files[c].write(0);
			}
			files[c].endln();
		}
	}
}

} //end namespaceRecal GenotypeLikelihoods
