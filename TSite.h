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

#define maxQualToPrint 1000
#define maxQualToPrintNaturalScale 1E-100

//---------------------------------------------------------------
//TSite
//---------------------------------------------------------------
class TSite{
protected:
	short int numGenotypes = 10;
	std::vector<TBase*>::iterator baseIterator;

	void normalizeGenotypeLikelihoods(double* emissionProbabilitiesPhredScaled, uint8_t* normalizedGL, uint32_t & maxLL, const int nGenotypes);

public:
	std::vector<TBase*> bases;
	bool hasData;
	double emissionProbabilities[10];
	char referenceBase; //optional

	TSite(){
		hasData = false;
		referenceBase = 'N';
	};

	TSite(TSite* other):TSite(){stealFromOther(other);};
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
	int refDepth();
	void addToBaseFrequencies(TBaseFrequencies & frequencies);
	double makePhred(double epsilon){
		return makePhredByRef(epsilon);
	};
	double makePhredByRef(double & epsilon){
		if(epsilon < maxQualToPrintNaturalScale) return maxQualToPrint;
		return -10.0 * log10(epsilon);
	};
	void calcEmissionProbabilities(double* vec);
	void calcEmissionProbabilities();
	void calculateP_g(double* genotypeProbabilities, double* P_g);
	double calculateWeightedSumOfEmissionProbs(double* weights);
	std::string getBases();
	std::string getEmissionProbs();
	double calculateLogLikelihood(double* genotypeProbabilities);
	void countAlleles(int* alleleCounts);
	void countAllelesForImbalance(long**** siteImbalance);
	void countMates(int* mateCounts);
	void countFwdRev(int* frCounts);
	void contextStats(int** contextCounts, TQualityMap & qualMap);
	void printPileup(gz::ogzstream & out);
	void printPileupToScreen();

	//MLE Callers
	void calculateDiploidPhredScaledGenotypeLikelihoods(double* phredGLs);
	void calculateNormalizedGenotypeLikelihoodsAndQuality(TRandomGenerator & randomGenerator, double* emissionProbabilitiesPhredScaled, double & quality, double & maxGenotypeProb, int & MLGenotype);
	void findSecondMostLikelyGenotype(TRandomGenerator & randomGenerator, double* emissionProbabilitiesPhredScaled, TGenotypeMap & genoMap, int MLGenotype, std::string & genoSecond);
	void findSecondMostProbableGenotype(TRandomGenerator & randomGenerator, double* postProb, TGenotypeMap & genoMap, int MAPGenotype, std::string & genoSecond);
	double calculatePHomozygous(double* pGenotype);

};

#endif /* TSITE_H_ */
