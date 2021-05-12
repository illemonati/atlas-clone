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
#include "TGenotypeMap.h"
#include "TReadGroups.h"

namespace GenotypeLikelihoods{

enum PMDTableType : uint8_t {forward3=0, forward5, reverse3, reverse5};

//---------------------------------------------------------------
//TPMDCounts
// counts are   [read group] [type] [from (ref)] [to (read)] [position]
// organized as TPMDTables -> TPMDTableReadGroup -> TPMDTPMDTable -> TPMDCounts -> countVec
//---------------------------------------------------------------

using countVec = std::vector<uint64_t>;

class TPMDCounts{
private:
	countVec _counts[4]; //_counts[A,C,G,T][position] are counts to (read)
	countVec _sums; //_sums[b] = sum_b _counts[b][p]
	uint16_t _size; //pos-specific for 0, 1, ..., _size - 2, then lumped into an extra bin (_size-1)
	uint16_t _sizeMinusOne;

	void _add(const uint16_t & pos, const BAM::Base & read);
	void _writeNormalizedOne(TOutputFile & out, countVec & these);

public:
	TPMDCounts(){
		_size = 0;
		_sizeMinusOne = 0;
	};
	TPMDCounts(const TPMDCounts & other) = default;
	~TPMDCounts() = default;

	uint16_t size() const { return _size; }
	void resize(const uint16_t & Size);
	void empty();
	void add(const uint16_t & pos, const BAM::Base & read);
	void add(const TPMDCounts & other);

	const countVec& operator[](const BAM::Base & b) const{
		return _counts[b];
	};
	const countVec& sums() const{
		return _sums;
	};

	void write(TOutputFile & out, const std::vector<std::string> & prefix, const bool & normalized);
};

//------------------------------------------------
// TPMDTable
//------------------------------------------------
class TPMDTable{
private:
	TPMDCounts _counts[4]; //_counts[A,C,G,T] are counts from (ref)

public:
	TPMDTable() = default;
	TPMDTable(const uint16_t & Size);
	TPMDTable(const TPMDTable & other);
	~TPMDTable() = default;

	uint16_t size() const { return _counts[0].size(); };
	void resize(const uint16_t & Size);
	void empty();
	void add(const uint16_t & pos, const BAM::Base & ref, const BAM::Base & read);
	void add(const TPMDTable & other);

	const TPMDCounts& operator[](const BAM::Base & b) const{
		return _counts[b];
	};

	void write(TOutputFile & out, std::vector<std::string> & prefix, const bool & normalized);
};

//------------------------------------------------
// TPMDTableReadGroup
//------------------------------------------------
class TPMDTableReadGroup{
private:
	TPMDTable _tables[4];

public:
	TPMDTableReadGroup(const uint16_t & TableLength);

	void add(const BAM::TSequencedBase & base, const BAM::Base & reference);

	const TPMDTable& operator[](const PMDTableType & Type) const{
		return _tables[Type];
	};

	void write(TOutputFile & out, std::vector<std::string> & prefix, const bool & normalized);
};

//------------------------------------------------
// TPMDTables
//------------------------------------------------
class TPMDTables{
private:
	const BAM::TReadGroupMap* _readGroupMap;
	const BAM::TReadGroups* _readGroups;
	uint16_t _tableLength;
	std::vector<TPMDTableReadGroup> _tables;

public:
	TPMDTables();
	TPMDTables(const BAM::TReadGroups* ReadGroups, const uint16_t & TableLength, const BAM::TReadGroupMap* ReadGroupMapObject);
	~TPMDTables(){};

	void initialize(const BAM::TReadGroups* ReadGroups, const uint16_t & TableLength, const BAM::TReadGroupMap* ReadGroupMap);
	const TPMDTableReadGroup& operator[](const uint16_t & ReadGroupID) const;

	void add(const BAM::TSequencedBase & base, const TGenotypeMap _genoMap;Base & reference);
	void write(std::string filename, const bool & normalize);
};

}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TPMDTABLES_H_ */
