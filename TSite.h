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
	short int numGenotypes;
	std::vector<TBase*>::iterator baseIterator;

public:
	std::vector<TBase*> bases;
	bool hasData;
	double* emissionProbabilities;
	char referenceBase; //optional

	TSite(){
		hasData = false;
		numGenotypes = 0;
		emissionProbabilities = NULL;
		referenceBase = 'N';
		//std::cout << "maxQualToPrintNaturalScale " << maxQualToPrintNaturalScale << std::endl;
	};
	TSite(TSite* other):TSite(){stealFromOther(other);};
	virtual ~TSite(){ clear(); };
	void clear();
	void stealFromOther(TSite* other);

	virtual void add(Base & base, int & quality, int PosInRead, int PosInReadRev, double thisPMD_CT, double thisPMD_GA, BaseContext & Context, int & ReadGroup){throw "Function 'add' Not implemented for base class TSite!"; };
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
	int depth(){return bases.size();};
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

	//MLE Callers
	void calculateNormalizedGenotypeLikelihoods(uint8_t* normalizedGL, uint32_t & maxLL);
	void calculateNormalizedGenotypeLikelihoodsAndQuality(TRandomGenerator & randomGenerator, double* emissionProbabilitiesPhredScaled, double & quality, double & maxGenotypeProb, int & MLGenotype);
	void findSecondMostLikelyGenotype(TRandomGenerator & randomGenerator, double* emissionProbabilitiesPhredScaled, TGenotypeMap & genoMap, int MLGenotype, std::string & genoSecond);
	void callMLEGenotype(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out);
	void callMLEGenotypeVCF(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool gVCF, bool noAltIfHomoRef, std::string & basesString);
	virtual void callMLEGenotypeKnownAlleles(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){ throw "callMLEGenotypeKnownAlleles not implemented for TSite base class!";};
	virtual void callMLEGenotypeVCFKnownAlleles(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt, bool noAltIfHomoRef, std::string & basesString){ throw "callMLEGenotypeVCFKnownAlleles not implemented for TSite base class!";};
	virtual void callMLEGenotypeKnownAllelesBeagle(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt, std::string & chr, int & pos, long & start, bool & printOnlyGL){ throw "callMLEGenotypeKnownAllelesBeagle not implemented for TSite base class!";};
	//Bayesian Callers
	void calculateGenotypePosteriorProbabilities(double* pGenotype, TRandomGenerator & randomGenerator, double* postProb, int & MAP);
	void callBayesianGenotype(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out);
	void callBayesianGenotypeVCF(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out);
	virtual void callBayesianGenotypeKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){ throw "callBayesianGenotypeKnownAlleles not implemented for TSite base class!";};
	virtual void callBayesianGenotypeVCFKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){ throw "callBayesianGenotypeVCFKnownAlleles not implemented for TSite base class!";};
	//Allele Presence Callers
	virtual void callAllelePresence(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out){ throw "callAllelePresence not implemented for TSite base class!";};
	virtual void callAllelePresenceKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){ throw "callAllelePresenceKnownAlleles not implemented for TSite base class!";};
	virtual void callAllelePresenceVCF(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool noAltIfHomoRef, std::string basesString){ throw "callAllelePresenceVCF not implemented for TSite base class!";};
	virtual void callAllelePresenceVCFKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt, bool noAltIfHomoRef, std::string basesString){ throw "callAllelePresenceVCFKnownAlleles not implemented for TSite base class!";};
	virtual void callRandomBase(TRandomGenerator & randomGenerator, gz::ogzstream & out) { throw "callRandomBase not implemented for TSite base class!";};
	virtual void majorityCall(TRandomGenerator & randomGenerator, gz::ogzstream & out){ throw "majorityCall not implemented for TSite base class!";};

	virtual double calculatePHomozygous(double* pGenotype){ throw "calculatePHomozygous not implemented for TSite base class!";};

	virtual void calculatePoolFreqLikelihoods(int & numChromosomes, TGenotypeMap & genoMap, Base & allele1, Base & allele2, gz::ogzstream & out){throw "calculatePoolFreqLikelihoods not implemented for TSite base class!";};
	void addToExpectedBaseCounts(TBaseFrequencies & baseFreq, double* expectedCounts){ throw "addToExpectedBaseCounts not implemented for TSite base class!";};
};

class TSiteDiploid:public TSite{
public:

	TSiteDiploid(){
		numGenotypes = 10;
		emissionProbabilities = new double[numGenotypes];
	};
	TSiteDiploid(TSite* other):TSiteDiploid(){stealFromOther(other);};
	~TSiteDiploid(){
		delete[] emissionProbabilities;
	};
	void add(Base & base, int & quality, int PosInRead, int PosInReadRev, double thisPMD_CT, double thisPMD_GA, BaseContext & Context, int & ReadGroup);

	//MLE callers
	void calculatePhredScaledGenotypeLikelihoodsKnownAlleles(TGenotypeMap & genoMap, char & alt, TRandomGenerator & randomGenerator, double* phredEmissionProbs, double & quality, double & maxGenotypeProb, int & MLGenotype);
	void calculateGenotypeLikelihoodsKnownAlleles(TGenotypeMap & genoMap, char & alt, TRandomGenerator & randomGenerator, double* emissionProbs, double & sumEmissionProbs, int & pos);
	void callMLEGenotypeKnownAlleles(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
	void callMLEGenotypeVCFKnownAlleles(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt, bool noAltIfHomoRef, std::string & basesString);
	void callMLEGenotypeKnownAllelesBeagle(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt, std::string & chr, int & pos, long & start, bool & printOnlyGL);
	//Bayesian Callers
	void calculateGenotypePosteriorProbabilitiesKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, char & alt, TRandomGenerator & randomGenerator, double* postProb, int & MAP);
	void callBayesianGenotypeKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
	void callBayesianGenotypeVCFKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
	//Allele Presence Callers
	void calculatePosteriorOnAllelePresence(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, double* postProbAllele, int & MAP);
	void callAllelePresence(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out);
	void callAllelePresenceVCF(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool noAltIfHomoRef, std::string basesString);
	void calculatePosteriorOnAllelePresenceKnownAlleles(double* pGenotype, char & alt, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, double* postProbAllele, int & MAP);
	void callAllelePresenceKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
	void callAllelePresenceVCFKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt, bool noAltIfHomoRef, std::string basesString);
	void callRandomBase(TRandomGenerator & randomGenerator, gz::ogzstream & out);
	void majorityCall(TRandomGenerator & randomGenerator, gz::ogzstream & out);
	double calculatePHomozygous(double* pGenotype);
};

class TSiteHaploid:public TSite{
public:

	TSiteHaploid(){
		numGenotypes = 4;
		emissionProbabilities = new double[numGenotypes];
	}
	TSiteHaploid(TSite* other):TSiteHaploid(){stealFromOther(other);};
	~TSiteHaploid(){
		delete[] emissionProbabilities;
	};
	void add(Base & base, int & quality, int PosInRead, int PosInReadRev, double thisPMD_CT, double thisPMD_GA, BaseContext & Context, int & ReadGroup);
	void addToExpectedBaseCounts(TBaseFrequencies & baseFreq, double* expectedCounts);
	void calculatePoolFreqLikelihoods(int & numChromosomes, TGenotypeMap & genoMap, Base & allele1, Base & allele2, gz::ogzstream & out);
};


#endif /* TSITE_H_ */
