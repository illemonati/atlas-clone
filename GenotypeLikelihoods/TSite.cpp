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

void TSite::add(TBase* base){
	bases.push_back(base);
	hasData = true;
};

void TSite::addToBaseFrequencies(TBaseFrequencies & frequencies){
	if(hasData){
		static double weight = 1.0 / bases.size();
		for(auto& b : bases){
			frequencies.add(b->data.base, weight);
		}
	}
};

std::string TSite::getBases(){
	if(!hasData) return "-";
	std::string b = "";
	for(auto& b : bases){
		b +=  getBaseAsChar(b->getBaseAsEnum());
	}
	return b;
};

unsigned int TSite::depth(){
	if(!hasData) return 0;
	return bases.size();
};

unsigned int TSite::refDepth(){
	if(!hasData) return 0;
	if(referenceBase == 'N') return 0;
	unsigned int counter = 0;
	for(unsigned int i=0; i<bases.size(); ++i){
		if(getBaseAsChar(bases[i]->getBaseAsEnum()) == referenceBase) ++counter;
	}
	return counter;
};

//PROBLEM
void TSite::calculateP_g(double* genotypeProbabilities, double* P_g){
	//calculate normalized genotype probabilities according to Bayes rule
	double sum = 0.0;
	for(int i=0; i<10; ++i){
		//P_g[i] =  emissionProbabilities[i] * genotypeProbabilities[i];
		sum += P_g[i];
	}
	for(int i=0; i<10; ++i){
		P_g[i] /= sum;
	}
};

void TSite::countAlleles(int* alleleCounts) const{
	alleleCounts[0] = 0;
	alleleCounts[1] = 0;
	alleleCounts[2] = 0;
	alleleCounts[3] = 0;

	for(TBase* it : bases)
		++alleleCounts[it->getBaseAsEnum()];
};

void TSite::countAllelesForImbalance(TAllelicDepthCounts & counts){
	//calculate and return imbalance
	int b[4] = {0};
	for(TBase* it : bases){
		++b[it->getBaseAsEnum()];
	}
	counts.addSite(b[0], b[1], b[2], b[3]);
};

void TSite::countMates(int* mateCounts){
	mateCounts[0] = 0;
	mateCounts[1] = 0;

	for(TBase* it : bases)
		++mateCounts[it->isSecondMate()];
};

void TSite::countFwdRev(int* frCounts){
	frCounts[0] = 0;
	frCounts[1] = 0;

	for(TBase* it : bases)
		++frCounts[it->isReverseStrand()];
};

void TSite::contextStats(int** contextCounts, TQualityMap & qualMap){
	for(TBase* it : bases){
		int q = qualMap.errorToPhredInt(it->errorRate);
		int c = it->data.context;
		++contextCounts[q][c];
	}
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


























