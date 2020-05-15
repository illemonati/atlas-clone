/*
 * TBase.h
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#ifndef TSITE_H_
#define TSITE_H_

#include <math.h>
#include <vector>
#include "TParameters.h"
#include "TGenotypeMap.h"
#include "TQualityMap.h"
#include "TBase.h"
#include "gzstream.h"
#include <algorithm>
#include "TRandomGenerator.h"
#include "../TAllelicDepthCounts.h"
#include "auxiliaryTools.h"

#define maxQualToPrint 1000
#define maxQualToPrintNaturalScale 1E-100

//---------------------------------------------------------------
//TSite
//---------------------------------------------------------------
class TSite{
protected:
	void normalizeGenotypeLikelihoods(double* emissionProbabilitiesPhredScaled, uint8_t* normalizedGL, uint32_t & maxLL, const int nGenotypes);

public:
	std::vector<TBase*> bases;
	bool hasData;
	char referenceBase; //optional

	GenotypeLikelihoods::TGenotypeLikelihoods genotypeLikelihoods;


	TSite(){
		hasData = false;
		referenceBase = 'N';
	};

	//TSite(TSite* other):TSite(){stealFromOther(other);};
	virtual ~TSite(){ clear(); };
	void clear();
	void stealFromOther(TSite* other);

	void add(TBase* base);
	void setRefBase(char & Base){
		if(Base == 'A' || Base == 'C' || Base == 'G' || Base == 'T' || Base == 'a' || Base == 'c' || Base == 'g' || Base == 't')
			referenceBase = Base;
		else referenceBase = 'N';
	};
	char getRefBase(){return referenceBase;};
	Base getRefBaseAsEnum(){
		if(referenceBase == 'A' || referenceBase == 'a') return A;
		if(referenceBase == 'C' || referenceBase == 'c') return C;
		if(referenceBase == 'G' || referenceBase == 'g') return G;
		if(referenceBase == 'T' || referenceBase == 't') return T;
		return N;
	};
	char getBaseAsChar(Base base){
		if(base == A) return 'A';
		if(base == C) return 'C';
		if(base == G) return 'G';
		if(base == T) return 'T';
		return 'N';

	}
	unsigned int depth();
	unsigned int refDepth();
	void addToBaseFrequencies(TBaseFrequencies & frequencies);
	void calculateP_g(double* genotypeProbabilities, double* P_g);
	std::string getBases();
	std::string getEmissionProbs();

	void countAlleles(int* alleleCounts) const;
	void countAllelesForImbalance(TAllelicDepthCounts & counts);
	void countMates(int* mateCounts);
	void countFwdRev(int* frCounts);
	void contextStats(int** contextCounts, TQualityMap & qualMap);
	void printPileup(TOutputFileZipped & out);

	//MLE Callers
	void calculateDiploidPhredScaledGenotypeLikelihoods(double* phredGLs);
	void calculateNormalizedGenotypeLikelihoodsAndQuality(TRandomGenerator & randomGenerator, double* emissionProbabilitiesPhredScaled, double & quality, double & maxGenotypeProb, int & MLGenotype);
	void findSecondMostLikelyGenotype(TRandomGenerator & randomGenerator, double* emissionProbabilitiesPhredScaled, TGenotypeMap & genoMap, int MLGenotype, std::string & genoSecond);
	void findSecondMostProbableGenotype(TRandomGenerator & randomGenerator, double* postProb, TGenotypeMap & genoMap, int MAPGenotype, std::string & genoSecond);
	double calculatePHomozygous(double* pGenotype);

};

#endif /* TSITE_H_ */
