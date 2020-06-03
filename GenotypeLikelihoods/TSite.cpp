/*
 * TBase.cpp
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#include "TSite.h"

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

void TSite::add(const TBase * base){
	bases.emplace_back(base);
	hasData = true;
};

void TSite::addToBaseFrequencies(TBaseFrequencies & frequencies){
	if(hasData){
		static double weight = 1.0 / bases.size();
		for(auto& b : bases){
			frequencies.add(b->base, weight);
		}
	}
};

std::string TSite::getBases(TGenotypeMap & genoMap){
	if(!hasData) return "-";
	std::string s = "";
	for(auto& b : bases){
		s +=  genoMap.getBaseAsChar(b->base);
	}
	return s;
};

uint32_t TSite::depth(){
	if(!hasData) return 0;
	return bases.size();
};

uint32_t TSite::refDepth(){
	if(!hasData) return 0;
	if(referenceBase == 'N') return 0;
	uint32_t counter = 0;
	for(auto& b : bases){
		if(b->base == referenceBase)
			++counter;
	}
	return counter;
};

void TSite::countAlleles(GenotypeLikelihoods::TBaseData & alleleCounts) const{
	alleleCounts[0] = 0;
	alleleCounts[1] = 0;
	alleleCounts[2] = 0;
	alleleCounts[3] = 0;

	for(auto& b : bases){
		++alleleCounts[b->base];
};

void TSite::countMates(int* mateCounts){
	mateCounts[0] = 0;
	mateCounts[1] = 0;

	for(TBaseOld* it : bases)
		++mateCounts[it->isSecondMate()];
};

void TSite::countFwdRev(int* frCounts){
	frCounts[0] = 0;
	frCounts[1] = 0;

	for(TBaseOld* it : bases)
		++frCounts[it->isReverseStrand()];
};

void TSite::printPileup(TOutputFileZipped & out){
	out << referenceBase;
	out << depth() << refDepth();
	out << getBases();
};


//-----------------------------------------------------------------------
// Genotype likelihoods
//-----------------------------------------------------------------------

//PROBLEM
double TSite::calculatePHomozygous(double* pGenotype){
	//calculate posterior probability for each genotype
	double postProb[10];
	double tot = 0.0;

	for(int i=0; i<10; ++i){
		//postProb[i] = emissionProbabilities[i] * pGenotype[i];
		tot += postProb[i];
	}

	//make sum for all homozygous genotypes
	return (postProb[AA] + postProb[CC] + postProb[GG] + postProb[TT]) / tot;
};


























