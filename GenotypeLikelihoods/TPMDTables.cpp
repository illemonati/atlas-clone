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
void TPMDCounts::_add(const uint16_t & pos, const Base & read){
	++_counts[read][pos];
	++_sums[pos];
};
void TPMDCounts::resize(const uint16_t & Size){
	_size = Size;
	_counts[A].resize(_size, 0);
	_counts[G].resize(_size, 0);
	_counts[C].resize(_size, 0);
	_counts[T].resize(_size, 0);
	_sums.resize(_size, 0);
};

void TPMDCounts::empty(){
	std::fill(_counts[A].begin(), _counts[A].end(), 0);
	std::fill(_counts[C].begin(), _counts[C].end(), 0);
	std::fill(_counts[G].begin(), _counts[G].end(), 0);
	std::fill(_counts[T].begin(), _counts[T].end(), 0);
	std::fill(_sums.begin(), _sums.end(), 0);
};

void TPMDCounts::add(const uint16_t & pos, const Base & read){
	if(pos < _sizeMinusOne){
		_add(pos, read);
	} else {
		_add(_sizeMinusOne, read);
	}
};

void TPMDCounts::add(const TPMDCounts & other){
	if(_size != other._size){
		for(uint16_t i = 0; i< _size; ++i){
			_counts[A][i] += other._counts[A][i];
			_counts[C][i] += other._counts[C][i];
			_counts[G][i] += other._counts[G][i];
			_counts[T][i] += other._counts[T][i];
			_sums[i] += other._sums[i];
		}
	}
};

void TPMDCounts::_writeNormalizedOne(TOutputFile & out, countVec & these){
	for(uint16_t i = 0; i< _size; ++i){
		out << (double) these[i] / (double) _sums[i];
	}
	out << std::endl;
};

void TPMDCounts::write(TOutputFile & out, const std::vector<std::string> & prefix, const bool & normalized){
	if(normalized){
		out << prefix << "A"; _writeNormalizedOne(out, _counts[A]);
		out << prefix << "C"; _writeNormalizedOne(out, _counts[C]);
		out << prefix << "G"; _writeNormalizedOne(out, _counts[G]);
		out << prefix << "T"; _writeNormalizedOne(out, _counts[T]);
	} else {
		out << prefix << "A" << _counts[A] << std::endl;
		out << prefix << "C" << _counts[C] << std::endl;
		out << prefix << "G" << _counts[G] << std::endl;
		out << prefix << "T" << _counts[T] << std::endl;
	}
};

//---------------------------------------------------------------
//TPMDTable
//---------------------------------------------------------------
TPMDTable::TPMDTable(const uint16_t & Size){
	resize(Size);
};

TPMDTable::TPMDTable(const TPMDTable & other){
	_counts[A] = other._counts[A];
	_counts[C] = other._counts[C];
	_counts[G] = other._counts[G];
	_counts[T] = other._counts[T];
};

void TPMDTable::resize(const uint16_t & Size){
	_counts[A].resize(Size);
	_counts[C].resize(Size);
	_counts[G].resize(Size);
	_counts[T].resize(Size);
};

void TPMDTable::empty(){
	_counts[A].empty();
	_counts[C].empty();
	_counts[G].empty();
	_counts[T].empty();
};

void TPMDTable::add(const uint16_t & pos, const Base & ref, const Base & read){
	_counts[ref].add(pos, read);
};

void TPMDTable::add(const TPMDTable & other){
	_counts[A].add(other[A]);
	_counts[C].add(other[C]);
	_counts[G].add(other[G]);
	_counts[T].add(other[T]);
};

void TPMDTable::write(TOutputFile & out, std::vector<std::string> & prefix, const bool & normalized){
	//add ref base to prefix
	prefix[3] = "A";
	_counts[A].write(out, prefix, normalized);

	prefix[3] = "C";
	_counts[C].write(out, prefix, normalized);

	prefix[3] = "G";
	_counts[G].write(out, prefix, normalized);

	prefix[3] = "T";
	_counts[T].write(out, prefix, normalized);
};

//------------------------------------------------
// TPMDTableReadGroup
//------------------------------------------------
void TPMDTableReadGroup::resize(const uint16_t & MaxLength){
	_tables[forward3].resize(MaxLength);
	_tables[forward5].resize(MaxLength);
	_tables[reverse3].resize(MaxLength);
	_tables[reverse5].resize(MaxLength);
};

void TPMDTableReadGroup::add(const BAM::TBase & base, const Base & reference){
	if(base.isReverseStrand()){
		Base ref = _genoMap.baseToFlippedBase[reference];
		Base bas = _genoMap.baseToFlippedBase[base.base];
		_tables[reverse3].add(base.distFrom3Prime, ref, bas);
		_tables[reverse5].add(base.distFrom5Prime, ref, bas);
	} else {
		_tables[forward3].add(base.distFrom3Prime, reference, base.base);
		_tables[forward5].add(base.distFrom5Prime, reference, base.base);
	}
};

void TPMDTableReadGroup::write(TOutputFile & out, std::vector<std::string> & prefix, const bool & normalized){
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
	readGroups = nullptr;
	maxReadLength = 0;
	readGroupMap = nullptr;
	origNumReadGroups = 0;
	numReadGroups = 0;
	_initialized = false;
};

TPMDTables::TPMDTables(BAM::TReadGroups* ReadGroups, int maxLengthForInference, int MaxReadLength, BAM::TReadGroupMap* ReadGroupMap){
	initialize(ReadGroups, maxLengthForInference, MaxReadLength, ReadGroupMap);
};

void TPMDTables::initialize(BAM::TReadGroups* ReadGroups, int tableLength, int MaxReadLength, BAM::TReadGroupMap* ReadGroupMap){
	readGroups = ReadGroups;
	maxReadLength = MaxReadLength;
	_tableLength = tableLength;
	readGroupMap = ReadGroupMap;
	origNumReadGroups = readGroupMap->getOrigNumReadGroups();
	numReadGroups = readGroupMap->getNumReadGroups();

	_tables.resize(numReadGroups, TPMDTableReadGroup(_tableLength));
};

const TPMDTableReadGroup& TPMDTables::operator[](const uint16_t & ReadGroupID) const{
	return _tables[readGroupMap->getIndex(ReadGroupID)];
};

void TPMDTables::add(const BAM::TBase & base, const Base & reference){
	_tables[readGroupMap->getIndex(base.readGroupID)].add(base, reference);
};

void TPMDTables::write(std::string filename, const bool & normalize){
	//compile header
	std::vector<std::string> header = {"ReadGroup", "direction", "fromEnd", "referenceBase", "sequencedBase" };
	for(uint16_t p = 1; p <= _tableLength; ++p){
		header.push_back("position_" + toString(p));
	}
	header.push_back("position_>" + toString(_tableLength));

	//open file
	TOutputFile out(filename, header);

	//loop over all read groups
	std::vector<std::string> prefix(4);
	for(int i=0; i<origNumReadGroups; ++i){
		if(readGroups->readGroupInUse(i)){
			uint16_t index = readGroupMap->getIndex(i);
			prefix[0] = readGroups->getName(i);
			_tables[index].write(out, prefix, normalize);
		}
	}
	out.close();
};

}; //end namespace
