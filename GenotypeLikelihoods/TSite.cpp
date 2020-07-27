/*
 * TBase.cpp
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#include "TSite.h"

namespace GenotypeLikelihoods{

//-------------------------------------------------------
//TSite
//-------------------------------------------------------
void TSite::clear(){
	_bases.clear();
	_referenceBase = N;
};

void TSite::add(BAM::TBase* base){
	_bases.push_back(base);
};

void TSite::addToBaseFrequencies(TBaseData & frequencies) const{
	if(!empty()){
		static double weight = 1.0 / _bases.size();
		for(auto& b : _bases){
			frequencies.add(b->base, weight);
		}
	}
};

std::string TSite::getBases(const TGenotypeMap & genoMap) const{
	if(empty()) return "-";
	std::string s = "";
	for(auto& b : _bases){
		s +=  genoMap.getBaseAsChar(b->base);
	}
	return s;
};

std::string TSite::getQualities(const TQualityMap & qualMap) const{
	if(empty()) return "-";
	std::string s = "";
	for(auto& b : _bases){
		s +=  (char) qualMap.phredIntToQuality(b->recalibratedQualityAsPhredInt);
	}
	return s;
};


uint32_t TSite::depth() const{
	return _bases.size();
};

uint32_t TSite::refDepth() const{
	if(empty()) return 0;
	if(_referenceBase == 'N') return 0;
	uint32_t counter = 0;
	for(auto& b : _bases){
		if(b->base == _referenceBase)
			++counter;
	}
	return counter;
};

void TSite::countAlleles(TBaseCounts & alleleCounts) const{
	alleleCounts.reset();
	for(auto& b : _bases){
		alleleCounts.add(b->base);
	}
};

void TSite::countMates(int* mateCounts) const{
	mateCounts[0] = 0;
	mateCounts[1] = 0;

	for(auto& it : _bases){
		++mateCounts[it->isSecondMate()];
	}
};

void TSite::countFwdRev(int* frCounts) const{
	frCounts[0] = 0;
	frCounts[1] = 0;

	for(auto& it : _bases){
		++frCounts[it->isReverseStrand()];
	}
};

//-------------------------------------------------------
//TSiteStorage
//-------------------------------------------------------
TSiteStorage::TSiteStorage(const TSite & site){
	_referenceBase = site.refBase();
	_genotype = site.genotype();
	for(std::vector<BAM::TBase*>::const_iterator b = site.cbegin(); b != site.cend(); ++b){
		_bases.emplace_back(**b);
	}
};

void TSiteStorage::clear(){
	_bases.clear();
	_referenceBase = N;
};

void TSiteStorage::add(const BAM::TBase * base){
	_bases.emplace_back(*base);
};

void TSiteStorage::add(const BAM::TBase & base){
	_bases.emplace_back(base);
};

void TSiteStorage::addToBaseFrequencies(TBaseData & frequencies) const{
	if(!empty()){
		static double weight = 1.0 / _bases.size();
		for(auto& b : _bases){
			frequencies.add(b.base, weight);
		}
	}
};

uint32_t TSiteStorage::depth() const{
	return _bases.size();
};


}; //end namespace






















