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
using countVec = std::vector<uint64_t>;
using PMDCounts = std::array<countVec, 4>;

//------------------------------------------------
// TPMDTable
//------------------------------------------------
class TPMDTable {
private:
	//TPMDCounts _counts[4]; 
	std::array<PMDCounts, 4> _counts; //_counts[A,C,G,T] are counts from (ref)
	std::array<countVec, 4> _sums;

public:
	TPMDTable() = default;
	TPMDTable(size_t Size);
	TPMDTable(const TPMDTable &other) = default;
	~TPMDTable() = default;

	size_t size() const { return _sums[0].size(); };
	void resize(size_t Size);
	void empty();
	void add(size_t pos, const genometools::Base &ref, const genometools::Base &read);
	void add(const TPMDTable &other);

	const PMDCounts &operator[](genometools::Base b) const { return _counts[b.get()]; };
	const countVec &sums(genometools::Base b) const { return _sums[b.get()]; }

	void write(coretools::TOutputFile &out, std::vector<std::string> &prefix, const bool &normalized);
};

//------------------------------------------------
// TPMDTableReadGroup
//------------------------------------------------
class TPMDTableReadGroup {
private:
	TPMDTable _tables[4];

public:
	TPMDTableReadGroup(uint16_t TableLength);

	void add(const BAM::TSequencedBase &base, const genometools::Base &reference);

	const TPMDTable &operator[](const PMDTableType &Type) const { return _tables[Type]; };

	void write(coretools::TOutputFile &out, std::vector<std::string> &prefix, const bool &normalized);
};

//------------------------------------------------
// TPMDTables
//------------------------------------------------
class TPMDTables {
private:
	const BAM::TReadGroupMap *_readGroupMap;
	const BAM::TReadGroups *_readGroups;
	uint16_t _tableLength;
	std::vector<TPMDTableReadGroup> _tables;

public:
	TPMDTables();
	TPMDTables(const BAM::TReadGroups *ReadGroups, uint16_t TableLength, const BAM::TReadGroupMap *ReadGroupMapObject);
	~TPMDTables(){};

	void initialize(const BAM::TReadGroups *ReadGroups, uint16_t TableLength, const BAM::TReadGroupMap *ReadGroupMap);
	const TPMDTableReadGroup &operator[](uint16_t ReadGroupID) const;

	void add(const BAM::TSequencedBase &base, const genometools::Base &reference);
	void write(std::string filename, const bool &normalize);
};

}; // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TPMDTABLES_H_ */
