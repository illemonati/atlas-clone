/*
 * TPMDTables.cpp
 *
 *  Created on: Apr 21, 2021
 *      Author: wegmannd
 */

#include "TPMDTables.h"

namespace GenotypeLikelihoods{

//---------------------------------------------------------------
//TPMDCounts
//---------------------------------------------------------------
void TPMDCounts::_add(const uint16_t & pos, const genometools::Base & read){
	++_counts[read.get()][pos];
	++_sums[pos];
};
void TPMDCounts::resize(const uint16_t & Size){
	_size = Size;
	_counts[BAM::A].resize(_size, 0);
	_counts[BAM::G].resize(_size, 0);
	_counts[BAM::C].resize(_size, 0);
	_counts[BAM::T].resize(_size, 0);
	_sums.resize(_size, 0);
};

void TPMDCounts::empty(){
	std::fill(_counts[BAM::A].begin(), _counts[BAM::A].end(), 0);
	std::fill(_counts[BAM::C].begin(), _counts[BAM::C].end(), 0);
	std::fill(_counts[BAM::G].begin(), _counts[BAM::G].end(), 0);
	std::fill(_counts[BAM::T].begin(), _counts[BAM::T].end(), 0);
	std::fill(_sums.begin(), _sums.end(), 0);
};

void TPMDCounts::add(const uint16_t & pos, const genometools::Base & read){
	if(pos < _sizeMinusOne){
		_add(pos, read);
	} else {
		_add(_sizeMinusOne, read);
	}
};

void TPMDCounts::add(const TPMDCounts & other){
	if(_size != other._size){
		for(uint16_t i = 0; i< _size; ++i){
			_counts[BAM::A][i] += other._counts[BAM::A][i];
			_counts[BAM::C][i] += other._counts[BAM::C][i];
			_counts[BAM::G][i] += other._counts[BAM::G][i];
			_counts[BAM::T][i] += other._counts[BAM::T][i];
			_sums[i] += other._sums[i];
		}
	}
};

void TPMDCounts::_writeNormalizedOne(coretools::TOutputFile & out, countVec & these){
	for(uint16_t i = 0; i< _size; ++i){
		out << (double) these[i] / (double) _sums[i];
	}
	out << std::endl;
};

void TPMDCounts::write(coretools::TOutputFile & out, const std::vector<std::string> & prefix, const bool & normalized){
	if(normalized){
		out << prefix << "A"; _writeNormalizedOne(out, _counts[BAM::A]);
		out << prefix << "C"; _writeNormalizedOne(out, _counts[BAM::C]);
		out << prefix << "G"; _writeNormalizedOne(out, _counts[BAM::G]);
		out << prefix << "T"; _writeNormalizedOne(out, _counts[BAM::T]);
	} else {
		out << prefix << "A" << _counts[BAM::A] << std::endl;
		out << prefix << "C" << _counts[BAM::C] << std::endl;
		out << prefix << "G" << _counts[BAM::G] << std::endl;
		out << prefix << "T" << _counts[BAM::T] << std::endl;
	}
};

//---------------------------------------------------------------
//TPMDTable
//---------------------------------------------------------------
TPMDTable::TPMDTable(const uint16_t & Size){
	resize(Size);
};

TPMDTable::TPMDTable(const TPMDTable & other){
	_counts[BAM::A] = other._counts[BAM::A];
	_counts[BAM::C] = other._counts[BAM::C];
	_counts[BAM::G] = other._counts[BAM::G];
	_counts[BAM::T] = other._counts[BAM::T];
};

void TPMDTable::resize(const uint16_t & Size){
	_counts[BAM::A].resize(Size);
	_counts[BAM::C].resize(Size);
	_counts[BAM::G].resize(Size);
	_counts[BAM::T].resize(Size);
};

void TPMDTable::empty(){
	_counts[BAM::A].empty();
	_counts[BAM::C].empty();
	_counts[BAM::G].empty();
	_counts[BAM::T].empty();
};

void TPMDTable::add(const uint16_t & pos, const genometools::Base & ref, const genometools::Base & read){
	_counts[ref.get()].add(pos, read);
};

void TPMDTable::add(const TPMDTable & other){
	_counts[BAM::A].add(other[BAM::A]);
	_counts[BAM::C].add(other[BAM::C]);
	_counts[BAM::G].add(other[BAM::G]);
	_counts[BAM::T].add(other[BAM::T]);
};

void TPMDTable::write(coretools::TOutputFile & out, std::vector<std::string> & prefix, const bool & normalized){
	//add ref base to prefix
	prefix[3] = "A";
	_counts[BAM::A].write(out, prefix, normalized);

	prefix[3] = "C";
	_counts[BAM::C].write(out, prefix, normalized);

	prefix[3] = "G";
	_counts[BAM::G].write(out, prefix, normalized);

	prefix[3] = "T";
	_counts[BAM::T].write(out, prefix, normalized);
};

//------------------------------------------------
// TPMDTableReadGroup
//------------------------------------------------
TPMDTableReadGroup::TPMDTableReadGroup(const uint16_t & TableLength){
	//resize
	_tables[forward3].resize(TableLength);
	_tables[forward5].resize(TableLength);
	_tables[reverse3].resize(TableLength);
	_tables[reverse5].resize(TableLength);
};

void TPMDTableReadGroup::add(const BAM::TSequencedBase & base, const genometools::Base & reference){
	if(base.isReverseStrand()){
		_tables[reverse3].add(base.distFrom3Prime, reference.flipped(), base.base.flipped());
		_tables[reverse5].add(base.distFrom5Prime, reference.flipped(), base.base.flipped());
	} else {
		_tables[forward3].add(base.distFrom3Prime, reference, base.base);
		_tables[forward5].add(base.distFrom5Prime, reference, base.base);
	}
};

void TPMDTableReadGroup::write(coretools::TOutputFile & out, std::vector<std::string> & prefix, const bool & normalized){
	prefix[1] = "forward";
	prefix[2] = "5'";
	_tables[forward5].write(out, prefix, normalized);

	prefix[2] = "3'";
	_tables[forward3].write(out, prefix, normalized);

	prefix[1] = "reverse";
	prefix[2] = "5'";
	_tables[reverse5].write(out, prefix, normalized);

	prefix[2] = "3'";
	_tables[reverse3].write(out, prefix, normalized);
};

//---------------------------------------------------------------
//TPMDTables
//---------------------------------------------------------------
TPMDTables::TPMDTables(){
	_tableLength = 0;
	_readGroups = nullptr;
	_readGroupMap = nullptr;
};

TPMDTables::TPMDTables(const BAM::TReadGroups* ReadGroups, const uint16_t & TableLength, const BAM::TReadGroupMap* ReadGroupMap){
	_tableLength = 0;
	initialize(ReadGroups, TableLength, ReadGroupMap);
};

void TPMDTables::initialize(const BAM::TReadGroups* ReadGroups, const uint16_t & TableLength, const BAM::TReadGroupMap* ReadGroupMap){
	if(_tableLength > 0){
		_tables.clear();
	}

	_readGroups = ReadGroups;
	_readGroupMap = ReadGroupMap;

	_tableLength = TableLength;
	_tables.resize(_readGroupMap->numReadGroupsInUse(), TPMDTableReadGroup(_tableLength));
};

const TPMDTableReadGroup& TPMDTables::operator[](const uint16_t & ReadGroupID) const{
	return _tables[_readGroupMap->pooledIndex(ReadGroupID)];
};

void TPMDTables::add(const BAM::TSequencedBase & base, const genometools::Base & reference){
	_tables[_readGroupMap->pooledIndex(base.readGroupID)].add(base, reference);
};

void TPMDTables::write(std::string filename, const bool & normalize){
	//compile header
	std::vector<std::string> header = {"ReadGroup", "direction", "fromEnd", "referenceBase", "sequencedBase" };
	for(uint16_t p = 1; p <= _tableLength; ++p){
		header.push_back("position_" + coretools::toString(p));
	}
	header.push_back("position_>" + coretools::toString(_tableLength));

	//open file
	coretools::TOutputFile out(filename, header);

	//loop over all read groups
	std::vector<std::string> prefix(4);
	for(int i=0; i<_readGroups->size(); ++i){
		if(_readGroups->readGroupInUse(i)){
			uint16_t index = _readGroupMap->pooledIndex(i);
			prefix[0] = _readGroups->getName(i);
			_tables[index].write(out, prefix, normalize);
		}
	}
	out.close();
};

}; //end namespace
