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
	if(hasData){
		bases.clear();
		hasData = false;
		referenceBase = 'N';
	}
};

/*
void TSite::stealFromOther(TSite* other){
	//this function extracts all data from the other object and sets it to empty
	hasData = other->hasData;
	if(hasData){
		//copy data
		referenceBase = other->referenceBase;
		genotypeLikelihoods = other->genotypeLikelihoods;

		//copy pointers to bases, BUT NOT BASES
		for(std::vector<TBase*>::iterator it = other->bases.begin(); it!=other->bases.end(); ++it){
			bases.push_back(*it);
		}
		//remove pointers from other site
		other->bases.clear();
		other->hasData = false;
	}
};
*/

void TSite::add(const BAM::TBase * base){
	bases.emplace_back(base);
	hasData = true;
};

void TSite::addToBaseFrequencies(TBaseData & frequencies) const{
	if(hasData){
		static double weight = 1.0 / bases.size();
		for(auto& b : bases){
			frequencies.add(b->base, weight);
		}
	}
};

std::string TSite::getBases(const TGenotypeMap & genoMap) const{
	if(!hasData) return "-";
	std::string s = "";
	for(auto& b : bases){
		s +=  genoMap.getBaseAsChar(b->base);
	}
	return s;
};

std::string TSite::getQualities(const TQualityMap & qualMap) const{
	if(!hasData) return "-";
	std::string s = "";
	for(auto& b : bases){
		s +=  (char) qualMap.phredIntToQuality(b->recalibratedQualityAsPhredInt);
	}
	return s;
};

uint32_t TSite::depth() const{
	if(!hasData) return 0;
	return bases.size();
};

uint32_t TSite::refDepth() const{
	if(!hasData) return 0;
	if(referenceBase == 'N') return 0;
	uint32_t counter = 0;
	for(auto& b : bases){
		if(b->base == referenceBase)
			++counter;
	}
	return counter;
};

void TSite::countAlleles(TBaseCounts & alleleCounts) const{
	alleleCounts.reset();
	for(auto& b : bases){
		alleleCounts.add(b->base);
	}
};

void TSite::countMates(int* mateCounts) const{
	mateCounts[0] = 0;
	mateCounts[1] = 0;

	for(auto& it : bases){
		++mateCounts[it->isSecondMate()];
	}
};

void TSite::countFwdRev(int* frCounts) const{
	frCounts[0] = 0;
	frCounts[1] = 0;

	for(auto& it : bases){
		++frCounts[it->isReverseStrand()];
	}
};


}; //end namespace






















