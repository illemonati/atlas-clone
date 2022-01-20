/*
 * TPMDTables.cpp
 *
 *  Created on: Apr 21, 2021
 *      Author: wegmannd
 */

#include "TPMDTables.h"

namespace GenotypeLikelihoods {

//---------------------------------------------------------------
// TPMDTable
//---------------------------------------------------------------
void TPMDTable::resize(size_t Size) {
	for (auto &from : _counts)
		for (auto & to: from)
			to.resize(Size);
	for (auto &s: _sums) s.resize(Size);
}

void TPMDTable::empty() {
	for (auto &from : _counts)
		for (auto & to: from)
			std::fill(to.begin(), to.end(), 0);
	for (auto &s: _sums) std::fill(s.begin(), s.end(), 0);
}

void TPMDTable::add(size_t pos, genometools::Base ref, genometools::Base read) {
	const auto p = std::min(pos, size() - 1);
	++_counts[ref.get()][read.get()][p];
	++_sums[ref.get()][p];
}

void TPMDTable::add(const TPMDTable &other) {
	if (size() != other.size()) return;

	for (size_t f = 0; f < _counts.size(); ++f) {
		for (size_t i = 0; i < size(); ++i) {
			_sums[f][i] += other._sums[f][i];
			for (size_t t = 0; t < _counts[f].size(); ++t) { _counts[f][t][i] += other._counts[f][t][i]; }
		}
	}
}

void TPMDTable::write(coretools::TOutputFile &out, std::vector<std::string> &prefix, bool normalized) {
	using namespace genometools;
	// add ref base to prefix
	for (Base f = Base::min(); f < Base::max(); ++f) {
		prefix[3] = (std::string)f;
		for (Base t = Base::min(); t < Base::max(); ++t) {
			out << prefix << std::string(t);
			if (normalized) {
				for (uint16_t i = 0; i < size(); ++i)
					out << static_cast<double>(_counts[f.get()][t.get()][i]) / _sums[f.get()][i];

			} else {
				out << _counts[f.get()][t.get()];
			}
			out << '\n';
		}
	}
}

//---------------------------------------------------------------
// TPMDTables
//---------------------------------------------------------------
TPMDTables::TPMDTables() {
	_tableLength  = 0;
	_readGroups   = nullptr;
	_readGroupMap = nullptr;
}

TPMDTables::TPMDTables(const BAM::TReadGroups *ReadGroups, size_t TableLength,
		       const BAM::TReadGroupMap *ReadGroupMap) {
	_tableLength = 0;
	initialize(ReadGroups, TableLength, ReadGroupMap);
}

void TPMDTables::initialize(const BAM::TReadGroups *ReadGroups, size_t TableLength,
			    const BAM::TReadGroupMap *ReadGroupMap) {
	if (_tableLength > 0) { _tables.clear(); }

	_readGroups   = ReadGroups;
	_readGroupMap = ReadGroupMap;

	_tableLength = TableLength;
	_tables.resize(_readGroupMap->numReadGroupsInUse());
	for (auto & ts: _tables)
		for (auto & t: ts) t.resize(TableLength);
}

void TPMDTables::add(const BAM::TSequencedBase &base, genometools::Base reference) {
	const auto from3 = base.distFrom3Prime < base.distFrom5Prime;
	const auto i_rg  = _readGroupMap->pooledIndex(base.readGroupID);
	if (base.isReverseStrand()) {
		if (from3)
			_tables[i_rg][reverse3].add(base.distFrom3Prime, reference.flipped(), base.base.flipped());
		else
			_tables[i_rg][reverse5].add(base.distFrom5Prime, reference.flipped(), base.base.flipped());
	} else {
		if (from3)
			_tables[i_rg][forward3].add(base.distFrom3Prime, reference, base.base);
		else
			_tables[i_rg][forward5].add(base.distFrom5Prime, reference, base.base);
	}
}

void TPMDTables::write(std::string filename, bool normalize) {
	// compile header
	std::vector<std::string> header = {"ReadGroup", "direction", "fromEnd", "referenceBase", "sequencedBase"};
	for (uint16_t p = 1; p <= _tableLength; ++p) { header.push_back("position_" + coretools::toString(p)); }
	header.push_back("position_>" + coretools::toString(_tableLength));

	// open file
	coretools::TOutputFile out(filename, header);

	// loop over all read groups
	static const std::array<std::string, 4> prefix1 = {"forward", "reverse", "forward", "reverse"};
	static const std::array<std::string, 4> prefix2 = {"5'", "3'", "5'", "3'"};
	std::vector<std::string> prefix(4);
	for (int i = 0; i < _readGroups->size(); ++i) {
		if (_readGroups->readGroupInUse(i)) {
			uint16_t index = _readGroupMap->pooledIndex(i);
			prefix[0]      = _readGroups->getName(i);
			for (size_t j = 0; j < 4; ++j) {
				prefix[1] =prefix1[j];
				prefix[2] =prefix2[j];
				_tables[index][j].write(out, prefix, normalize);
			}
		}
	}
	out.close();
};

}; // namespace GenotypeLikelihoods
