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
	_counts[genometools::A].resize(_size, 0);
	_counts[genometools::G].resize(_size, 0);
	_counts[genometools::C].resize(_size, 0);
	_counts[genometools::T].resize(_size, 0);
	_sums.resize(_size, 0);
};

void TPMDCounts::empty(){
	std::fill(_counts[genometools::A].begin(), _counts[genometools::A].end(), 0);
	std::fill(_counts[genometools::C].begin(), _counts[genometools::C].end(), 0);
	std::fill(_counts[genometools::G].begin(), _counts[genometools::G].end(), 0);
	std::fill(_counts[genometools::T].begin(), _counts[genometools::T].end(), 0);
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
			_counts[genometools::A][i] += other._counts[genometools::A][i];
			_counts[genometools::C][i] += other._counts[genometools::C][i];
			_counts[genometools::G][i] += other._counts[genometools::G][i];
			_counts[genometools::T][i] += other._counts[genometools::T][i];
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
		out << prefix << "A"; _writeNormalizedOne(out, _counts[genometools::A]);
		out << prefix << "C"; _writeNormalizedOne(out, _counts[genometools::C]);
		out << prefix << "G"; _writeNormalizedOne(out, _counts[genometools::G]);
		out << prefix << "T"; _writeNormalizedOne(out, _counts[genometools::T]);
	} else {
		out << prefix << "A" << _counts[genometools::A] << std::endl;
		out << prefix << "C" << _counts[genometools::C] << std::endl;
		out << prefix << "G" << _counts[genometools::G] << std::endl;
		out << prefix << "T" << _counts[genometools::T] << std::endl;
	}
};

//---------------------------------------------------------------
//TPMDTable
//---------------------------------------------------------------
TPMDTable::TPMDTable(const uint16_t & Size){
	resize(Size);
};

TPMDTable::TPMDTable(const TPMDTable & other){
	_counts[genometools::A] = other._counts[genometools::A];
	_counts[genometools::C] = other._counts[genometools::C];
	_counts[genometools::G] = other._counts[genometools::G];
	_counts[genometools::T] = other._counts[genometools::T];
};

void TPMDTable::resize(const uint16_t & Size){
	_counts[genometools::A].resize(Size);
	_counts[genometools::C].resize(Size);
	_counts[genometools::G].resize(Size);
	_counts[genometools::T].resize(Size);
};

void TPMDTable::empty(){
	_counts[genometools::A].empty();
	_counts[genometools::C].empty();
	_counts[genometools::G].empty();
	_counts[genometools::T].empty();
};

void TPMDTable::add(const uint16_t & pos, const genometools::Base & ref, const genometools::Base & read){
	_counts[ref.get()].add(pos, read);
};

void TPMDTable::add(const TPMDTable & other){
	_counts[genometools::A].add(other[genometools::A]);
	_counts[genometools::C].add(other[genometools::C]);
	_counts[genometools::G].add(other[genometools::G]);
	_counts[genometools::T].add(other[genometools::T]);
};

void TPMDTable::write(coretools::TOutputFile & out, std::vector<std::string> & prefix, const bool & normalized){
	//add ref base to prefix
	prefix[3] = "A";
	_counts[genometools::A].write(out, prefix, normalized);

	prefix[3] = "C";
	_counts[genometools::C].write(out, prefix, normalized);

	prefix[3] = "G";
	_counts[genometools::G].write(out, prefix, normalized);

	prefix[3] = "T";
	_counts[genometools::T].write(out, prefix, normalized);
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
	//  if distFrom3 < distFrom 5
		_tables[reverse3].add(base.distFrom3Prime, reference.flipped(), base.base.flipped());
	//  else
		_tables[reverse5].add(base.distFrom5Prime, reference.flipped(), base.base.flipped());
	} else {
	//  if distFrom3 < distFrom 5
		_tables[forward3].add(base.distFrom3Prime, reference, base.base);
	//  else
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
