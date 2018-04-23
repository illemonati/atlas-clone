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
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>
#include "TParameters.h"
#include "TGenotypeMap.h"
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

//	void add(Base & base, int & quality, int PosInRead, int PosInReadRev, double thisPMD_CT, double thisPMD_GA, BaseContext & Context, int & ReadGroup);
	void add(TBase* base);
	void setRefBase(char & Base){
		if(Base == 'A' || Base == 'C' || Base == 'G' || Base == 'T')
			referenceBase = Base;
		else referenceBase = 'N';
	};
	char getRefBase(){return referenceBase;};
	Base getRefBaseAsEnum(){
		if(referenceBase == 'A') return A;
		if(referenceBase == 'C') return C;
		if(referenceBase == 'G') return G;
		if(referenceBase == 'T') return T;
		return N;
	};
	char getBaseAsChar(Base base){
		if(base == A) return 'A';
		if(base == C) return 'C';
		if(base == G) return 'G';
		if(base == T) return 'T';
		return 'N';

	}
	int depth();
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
	void calculateP_g(double* & genotypeProbabilities, double* & P_g);
	double calculateWeightedSumOfEmissionProbs(double* weights);
	std::string getBases();
	std::string getEmissionProbs();
	double calculateLogLikelihood(double* genotypeProbabilities);
	void countAlleles(long**** siteImbalance);
	void printPileup(gz::ogzstream & out);

	//MLE Callers
	void calculateNormalizedGenotypeLikelihoods(uint8_t* normalizedGL, uint32_t & maxLL);
	void calculateNormalizedGenotypeLikelihoodsAndQuality(TRandomGenerator & randomGenerator, double* emissionProbabilitiesPhredScaled, double & quality, double & maxGenotypeProb, int & MLGenotype);
	void findSecondMostLikelyGenotype(TRandomGenerator & randomGenerator, double* emissionProbabilitiesPhredScaled, TGenotypeMap & genoMap, int MLGenotype, std::string & genoSecond);
	void callMLEGenotype(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out);
	void callMLEGenotypeVCF(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool gVCF, bool noAltIfHomoRef, std::string & basesString);
	void callMLEGenotypeKnownAlleles(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
	void callMLEGenotypeVCFKnownAlleles(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt, bool noAltIfHomoRef, std::string & basesString);
	void calculateGenotypeLikelihoodsKnownAlleles(TGenotypeMap & genoMap, char & alt, TRandomGenerator & randomGenerator, double* emissionProbs, double & sumEmissionProbs, int & pos);
	void calculatePhredScaledGenotypeLikelihoodsKnownAlleles(TGenotypeMap & genoMap, char & alt, TRandomGenerator & randomGenerator, double* phredEmissionProbs, double & quality, double & maxGenotypeProb, int & MLGenotype);
	void callMLEGenotypeKnownAllelesBeagle(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt, std::string & chr, int & pos, long & start, bool & printOnlyGL);
	//Bayesian Callers
	void calculateGenotypePosteriorProbabilities(double* pGenotype, TRandomGenerator & randomGenerator, double* postProb, int & MAP);
	void callBayesianGenotype(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out);
	void callBayesianGenotypeVCF(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out);
	void callBayesianGenotypeKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
	void calculateGenotypePosteriorProbabilitiesKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, char & alt, TRandomGenerator & randomGenerator, double* postProb, int & MAP);
	void callBayesianGenotypeVCFKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
	//Allele Presence Callers
	void callAllelePresence(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out);
	void callAllelePresenceKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
	void calculatePosteriorOnAllelePresence(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, double* postProbAllele, int & MAP);
	void calculatePosteriorOnAllelePresenceKnownAlleles(double* pGenotype, char & alt, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, double* postProbAllele, int & MAP);
	void callAllelePresenceVCF(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool noAltIfHomoRef, std::string basesString);
	void callAllelePresenceVCFKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt, bool noAltIfHomoRef, std::string basesString);
	void callRandomBase(TRandomGenerator & randomGenerator, gz::ogzstream & out);
	void majorityCall(TRandomGenerator & randomGenerator, gz::ogzstream & out);
	double calculatePHomozygous(double* pGenotype);

};

#endif /* TSITE_H_ */
