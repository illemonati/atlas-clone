/*
 * TPMDTables.h
 *
 *  Created on: Apr 21, 2021
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TPMDTABLES_H_
#define GENOTYPELIKELIHOODS_TPMDTABLES_H_

#include "TBase.h"
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
	uint16_t _size; //pos-specific up 0, 1, ..., _size - 1, then lumped into an extra bin

	void _add(const uint16_t & pos, const Base & read);
	void _writeNormalizedOne(TOutputFile & out, countVec & these);

public:
	TPMDCounts(){
		_size = 0;
	};
	~TPMDCounts() = default;

	uint16_t size() const { return _size; }
	void resize(const uint16_t & Size);
	void empty();
	void add(const uint16_t & pos, const Base & read);

	const countVec& operator[](const Base & b) const{
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
	TPMDTable();
	TPMDTable(const uint16_t & Size);
	~TPMDTable() = default;

	uint16_t size() const { return _counts[0].size(); };
	void resize(const uint16_t & Size);
	void empty();
	void add(const uint16_t & pos, const Base & ref, const Base & read);

	const TPMDCounts& operator[](const Base & b) const{
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
	TGenotypeMap _genoMap;

public:
	TPMDTableReadGroup() = default;

	void resize(const uint16_t & MaxLength);
	void empty();
	void add(const BAM::TBase & base, const Base & reference);

	void write(TOutputFile & out, std::vector<std::string> & prefix, const bool & normalized);
};

//------------------------------------------------
// TPMDTables
//------------------------------------------------
class TPMDTables{
private:
	BAM::TReadGroupMap* readGroupMap;
	BAM::TReadGroups* readGroups;
	int maxReadLength;
	uint16_t _tableLength;
	int origNumReadGroups;
	int numReadGroups;
	std::vector<TPMDTableReadGroup> _tables;
	bool _initialized;

public:
	TPMDTables();
	TPMDTables(BAM::TReadGroups* ReadGroups, int tableLength, int MaxReadLength, BAM::TReadGroupMap* ReadGroupMapObject);
	~TPMDTables(){};

	void initialize(BAM::TReadGroups* ReadGroups, int maxLengthForInference, int MaxReadLength, BAM::TReadGroupMap* ReadGroupMapObject);

	void add(const BAM::TBase & base, const Base & reference);
	void write(std::string filename, const bool & normalize);
};

}; //end namespace

#endif /* GENOTYPELIKELIHOODS_TPMDTABLES_H_ */
