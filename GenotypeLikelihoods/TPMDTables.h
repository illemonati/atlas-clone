/*
 * TPMDTables.h
 *
 *  Created on: Apr 21, 2021
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TPMDTABLES_H_
#define GENOTYPELIKELIHOODS_TPMDTABLES_H_

#include <stddef.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "genometools/GenotypeTypes.h"
#include "TReadGroups.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/devtools.h"

namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }
namespace coretools { class TOutputFile; }

namespace GenotypeLikelihoods {



//------------------------------------------------
// TPMDTable
//------------------------------------------------
class TPMDTable {
private:
	static constexpr size_t _N = coretools::index(genometools::Base::max) + 1;
	std::vector<coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base, _N>, genometools::Base, _N>> _counts;

public:
	TPMDTable() = default;

	TPMDTable(size_t Size) { resize(Size); };
	size_t size() const { return _counts.size(); };
	void resize(size_t Size);
	const auto &operator[](size_t i) const noexcept {return _counts[i];}
	auto &operator[](size_t i) noexcept {return _counts[i];}
	void add(size_t pos, genometools::Base ref, genometools::Base read);
	void write(coretools::TOutputFile &out, std::array<std::string, 4> &prefix, bool normalized);
};

enum class ReadEnd : size_t { min = 0, forward5 = min, reverse5, forward3, reverse3, max };
using PMDTable_RG = coretools::TStrongArray<TPMDTable, ReadEnd>;

//------------------------------------------------
// TPMDTables
//------------------------------------------------
class TPMDTables {
private:
	const BAM::TReadGroupMap *_readGroupMap = nullptr;
	const BAM::TReadGroups *_readGroups = nullptr;
	size_t _tableLength = 0;
	std::vector<PMDTable_RG> _tables;
public:
	TPMDTables() = default;

	TPMDTables(const BAM::TReadGroups *ReadGroups, size_t TableLength, const BAM::TReadGroupMap *ReadGroupMap) {
		initialize(ReadGroups, TableLength, ReadGroupMap);
	}
	void initialize(const BAM::TReadGroups *ReadGroups, size_t TableLength, const BAM::TReadGroupMap *ReadGroupMap);
	const PMDTable_RG &operator[](size_t ReadGroupID) const { return _tables[_readGroupMap->pooledIndex(ReadGroupID)]; }
	void add(const BAM::TSequencedBase &base, genometools::Base reference);
	void write(std::string filename, bool normalize);
};

} // namespace GenotypeLikelihoods
#endif /* GENOTYPELIKELIHOODS_TPMDTABLES_H_ */
