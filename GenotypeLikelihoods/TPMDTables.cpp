/*
 * TPMDTables.cpp
 *
 *  Created on: Apr 21, 2021
 *      Author: wegmannd
 */

#include "TPMDTables.h"

#include <stdint.h>
#include <ostream>

#include "coretools/Files/TOutputFile.h"
#include "TReadGroups.h"
#include "TSequencedBase.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenotypeLikelihoods {

//---------------------------------------------------------------
// TPMDTable
//---------------------------------------------------------------
void TPMDTable::resize(size_t Size) {
	for (auto &from : _counts)
		for (auto &to : from) to.resize(Size + 1); // Size + 1 for all larger than Size
	for (auto &s : _sums) s.resize(Size + 1);          // Size + 1 for all larger than Size
}

void TPMDTable::empty() {
	for (auto &from : _counts)
		for (auto & to: from)
			std::fill(to.begin(), to.end(), 0);
	for (auto &s: _sums) std::fill(s.begin(), s.end(), 0);
}

void TPMDTable::add(size_t pos, genometools::Base ref, genometools::Base read) {
	const auto p = std::min(pos, size());
	++_counts[ref][read][p];
	++_sums[ref][p];
}

void TPMDTable::add(const TPMDTable &other) {
	using genometools::Base;
	if (size() != other.size()) return;

	//for (size_t f = 0; f < _counts.size(); ++f)
	for (Base f = Base::min; f < Base::max; ++f)
		for (Base t = Base::min; t < Base::max; ++t)
			for (size_t i = 0; i < size(); ++i) _counts[f][t][i] += other._counts[f][t][i];

	for (Base f = Base::min; f < Base::max; ++f)
		for (size_t i = 0; i < _sums[f].size(); ++i) 
			_sums[f][i] += other._sums[f][i];
}

void TPMDTable::write(coretools::TOutputFile &out, std::vector<std::string> &prefix, bool normalized) {
	using namespace genometools;
	for (Base f = Base::min; f < Base::max; ++f) {
		prefix[3] = toString(f);
		for (Base t = Base::min; t < Base::max; ++t) {
			out.write(prefix, toString(t));
			if (normalized) {
				for (uint16_t i = 0; i < _sums[f].size(); ++i)
					out.write(static_cast<double>(_counts[f][t][i])/_sums[f][i]);
			} else {
				out.write(_counts[f][t]);
			}
			out.endln();
		}
	}
}

//---------------------------------------------------------------
// TPMDTables
//---------------------------------------------------------------

void TPMDTables::initialize(const BAM::TReadGroups *ReadGroups, size_t TableLength,
			    const BAM::TReadGroupMap *ReadGroupMap) {
	if (_tableLength > 0) _tables.clear();

	_readGroups   = ReadGroups;
	_readGroupMap = ReadGroupMap;
	_tableLength  = TableLength;
	//_tables.resize(_readGroupMap->numReadGroupsInUse()); // ReadgroupMap not continuous
	_tables.resize(_readGroups->size());
	for (auto i = _readGroups->cbegin(); i != _readGroups->cend(); ++i) {
		auto & table = _tables[_readGroupMap->pooledIndex(i->id())];
		for (auto &t : table) t.resize(TableLength);
	}
}

void TPMDTables::add(const BAM::TSequencedBase &base, genometools::Base reference) {
	const auto from3 = base.distFrom3Prime < base.distFrom5Prime;
	auto & table = _tables[_readGroupMap->pooledIndex(base.readGroupID)];
	if (base.isReverseStrand()) {
		if (from3)
			table[reverse3].add(base.distFrom3Prime, flipped(reference), flipped(base.base));
		else
			table[reverse5].add(base.distFrom5Prime, flipped(reference), flipped(base.base));
	} else {
		if (from3)
			table[forward3].add(base.distFrom3Prime, reference, base.base);
		else
			table[forward5].add(base.distFrom5Prime, reference, base.base);
	}
}

void TPMDTables::write(std::string filename, bool normalize) {
	// compile header
	std::vector<std::string> header = {"readGroup", "direction", "fromEnd", "referenceBase", "sequencedBase"};
	for (size_t p = 1; p <= _tableLength; ++p) header.push_back("position_" + coretools::str::toString(p));
	header.push_back("position_>" + coretools::str::toString(_tableLength));

	coretools::TOutputFile out(filename, header);
	static const std::array<std::string, 4> prefix1 = {"forward", "reverse", "forward", "reverse"};
	static const std::array<std::string, 4> prefix2 = {"5'", "3'", "5'", "3'"};

	// loop over all read groups
	std::vector<std::string> prefix(4);
	for (int i = 0; i < _readGroups->size(); ++i) {
		if (_readGroups->readGroupInUse(i)) {
			auto & table = _tables[_readGroupMap->pooledIndex(i)];
			prefix[0]   = _readGroups->getName(i);
			for (size_t j = 0; j < 4; ++j) {
				prefix[1] = prefix1[j];
				prefix[2] = prefix2[j];
				table[j].write(out, prefix, normalize);
			}
		}
	}
	out.close();
}

} // namespace GenotypeLikelihoods
