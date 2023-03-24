/*
 * TPMDTables.cpp
 *
 *  Created on: Apr 21, 2021
 *      Author: wegmannd
 */

#include "TPMDTables.h"

#include <stdint.h>
#include <ostream>

#include "coretools/Containers/TStrongArray.h"
#include "coretools/Files/TOutputFile.h"
#include "TReadGroups.h"
#include "TSequencedBase.h"
#include "coretools/Strings/stringFunctions.h"
#include "genometools/GenotypeTypes.h"

namespace GenotypeLikelihoods {

//---------------------------------------------------------------
// TPMDTable
//---------------------------------------------------------------
void TPMDTable::resize(size_t Size) {
	_counts.resize(Size + 1);
}

void TPMDTable::add(size_t pos, genometools::Base ref, genometools::Base read) {
	const auto p = std::min(pos, size() - 1);
	++_counts[p][ref][read];
}

void TPMDTable::write(coretools::TOutputFile &out, std::array<std::string, 4> &prefix, bool normalized) {
	using namespace genometools;
	for (Base f = Base::min; f <= Base::max; ++f) {
		prefix[3] = toString(f);
		std::vector<size_t> sums(size(), 0.);
		for (size_t i = 0; i < size(); ++i) {
			for (Base t = Base::min; t < Base::max; ++t) { sums[i] += _counts[i][f][t]; }
		}

		for (Base t = Base::min; t <= Base::max; ++t) {
			out.write(prefix, toString(t));
			for (size_t i = 0; i < size(); ++i) {
				if (normalized)
					out.write(static_cast<double>(_counts[i][f][t]) / sums[i]);
				else
					out.write(_counts[i][f][t]);
			}
			out.endln();
		}
		if (!normalized) { out.writeln(prefix, "sum", sums); }
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
			table[ReadEnd::reverse3].add(base.distFrom3Prime, flipped(reference), flipped(base.base));
		else
			table[ReadEnd::reverse5].add(base.distFrom5Prime, flipped(reference), flipped(base.base));
	} else {
		if (from3)
			table[ReadEnd::forward3].add(base.distFrom3Prime, reference, base.base);
		else
			table[ReadEnd::forward5].add(base.distFrom5Prime, reference, base.base);
	}
}

void TPMDTables::write(std::string filename, bool normalize) {
	// compile header
	std::vector<std::string> header = {"readGroup", "strand", "fromEnd", "ref", "data"};
	for (size_t p = 1; p <= _tableLength; ++p) header.push_back("pos_" + coretools::str::toString(p));
	header.push_back("pos>" + coretools::str::toString(_tableLength));

	coretools::TOutputFile out(filename, header);

	constexpr auto directions = [](){
		coretools::TStrongArray<std::array<std::string_view, 2>, ReadEnd> ar{};
		ar[ReadEnd::forward3] = {"forward", "3'"};
		ar[ReadEnd::forward5] = {"forward", "5'"};
		ar[ReadEnd::reverse3] = {"reverse", "3'"};
		ar[ReadEnd::reverse5] = {"reverse", "5'"};
		return ar;
	}();

	// loop over all read groups
	std::array<std::string, 4> prefix;
	for (size_t i = 0; i < _readGroups->size(); ++i) {
		if (_readGroups->readGroupInUse(i)) {
			auto & table = _tables[_readGroupMap->pooledIndex(i)];
			prefix[0]    = _readGroups->getName(i);
			for (auto j = ReadEnd::min; j < ReadEnd::max; ++j) {
				prefix[1] = directions[j].front();
				prefix[2] = directions[j].back();
				table[j].write(out, prefix, normalize);
			}
		}
	}
	out.close();
}

} // namespace GenotypeLikelihoods
