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
/*
TSite::TSite(const TSite & site){
    _referenceBase = site.refBase();
    _genotype = site.genotype();
    for(std::vector<BAM::TBase*>::const_iterator b = site.cbegin(); b != site.cend(); ++b){
        _bases.emplace_back(**b);
}
};*/

void TSite::clear(){
	_bases.clear();
	_referenceBase = BAM::N;
};

void TSite::add(const BAM::TSequencedBase & base){
    _bases.emplace_back(base);
};

void TSite::addToBaseFrequencies(TBaseData & frequencies) const{
	if(!empty()){
		static double weight = 1.0 / _bases.size();
		for(auto& b : _bases){
			frequencies.add(b.base, weight);
		}
	}
};

void TSite::downsample(const uint32_t & maxDepth, const TSubsamplePicker & picker){
	//only subsample if depth > maxDepth
	if(_bases.size() > maxDepth){
		//select subsample
		const auto& subsample = picker.pick(_bases.size(), maxDepth);

		//copy bases to new vector
		std::vector<BAM::TSequencedBase> newBases;
		for(auto& it : subsample){
			newBases.emplace_back(_bases[it]);
		}

		//swap vectors
		_bases = newBases;
	}
};

std::string TSite::getBases() const{
	if(empty()) return "-";
	std::string s = "";
	for(auto& b : _bases){
		s +=  (char) b.base;
	}
	return s;
};

std::string TSite::getQualities() const{
	if(empty()) return "-";
	std::string s = "";
	for(auto& b : _bases){
		s +=  (char) BAM::BaseQuality(b.recalibratedQualityAsPhredInt);
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
		if(b.base == _referenceBase)
			++counter;
	}
	return counter;
};

void TSite::countAlleles(TBaseCounts & alleleCounts) const{
	alleleCounts.reset();
	for(auto& b : _bases){
		alleleCounts.add(b.base);
	}
};

void TSite::countMates(int* mateCounts) const{
	mateCounts[0] = 0;
	mateCounts[1] = 0;

	for(auto& it : _bases){
		++mateCounts[it.isSecondMate()];
	}
};

void TSite::countFwdRev(int* frCounts) const{
	frCounts[0] = 0;
	frCounts[1] = 0;

	for(auto& it : _bases){
		++frCounts[it.isReverseStrand()];
	}
};

}; //end namespace






















