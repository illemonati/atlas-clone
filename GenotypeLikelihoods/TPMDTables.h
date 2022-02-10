/*
 * TPMDTables.h
 *
 *  Created on: Apr 21, 2021
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TPMDTABLES_H_
#define GENOTYPELIKELIHOODS_TPMDTABLES_H_

#include "../BAM/TSequencedBase.h"
#include "TFile.h"
#include "TReadGroups.h"

namespace GenotypeLikelihoods {

enum PMDTableType : uint8_t { forward3 = 0, forward5, reverse3, reverse5 };

using countVec    = std::vector<uint64_t>;
using PMDCounts   = std::array<countVec, 4>;

//------------------------------------------------
// TPMDTable
//------------------------------------------------
class TPMDTable {
private:
	std::array<PMDCounts, 4> _counts; //_counts[A,C,G,T] are counts from (ref)
	std::array<countVec, 4> _sums;
public:
	TPMDTable() = default;
	TPMDTable(const TPMDTable &) = default;
	TPMDTable& operator=(const TPMDTable &) = default;
	~TPMDTable() = default;

	TPMDTable(size_t Size) { resize(Size); };
	size_t size() const { return _sums[0].size(); };
	void resize(size_t Size);
	void empty();
	void add(size_t pos, genometools::Base ref, genometools::Base read);
	void add(const TPMDTable &other);
	const PMDCounts &operator[](genometools::Base b) const { return _counts[genometools::index(b)]; }
	const countVec &sums(genometools::Base b) const { return _sums[genometools::index(b)]; }
	void write(coretools::TOutputFile &out, std::vector<std::string> &prefix, bool normalized);
};

using PMDTable_RG = std::array<TPMDTable, 4>;

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
	TPMDTables(const TPMDTables &) = default;
	TPMDTables& operator=(const TPMDTables &) = default;
	~TPMDTables() = default;

	TPMDTables(const BAM::TReadGroups *ReadGroups, size_t TableLength, const BAM::TReadGroupMap *ReadGroupMap) {
		initialize(ReadGroups, TableLength, ReadGroupMap);
	}
	void initialize(const BAM::TReadGroups *ReadGroups, size_t TableLength, const BAM::TReadGroupMap *ReadGroupMap);
	const PMDTable_RG &operator[](size_t ReadGroupID) const { return _tables[ReadGroupID]; }
	void add(const BAM::TSequencedBase &base, genometools::Base reference);
	void write(std::string filename, bool normalize);
};

} // namespace GenotypeLikelihoods
#endif /* GENOTYPELIKELIHOODS_TPMDTABLES_H_ */
