/*
 * TBase.h
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#ifndef TSITE_H_
#define TSITE_H_

#include "stringFunctions.h"
#include <math.h>
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>
#include "TParameters.h"
#include "TGenotypeMap.h"
#include "TBase.h"
#include "gzstream.h"
#include <algorithm>
#include "TRandomGenerator.h"


//---------------------------------------------------------------
//TSite
//---------------------------------------------------------------
class TSite{
public:
	bool hasData;
	std::vector<TBase*> bases;
	int numGenotypes;
	double* emissionProbabilities;
	double* P_g; //P(g|d, theta, pi), see equation (3)
	char referenceBase; //optional
	double maxQualToPrint, maxQualToPrintNaturalScale;

	TSite(){
		hasData = false;
		numGenotypes = 0;
		emissionProbabilities = NULL;
		P_g = NULL;
		referenceBase = 'N';
		maxQualToPrint = 1000;
		maxQualToPrintNaturalScale = pow(10.0, -maxQualToPrint / 10.0);

	};
	virtual ~TSite(){ clear(); };

	void clear();

	virtual void add(char & base, char & quality, int PosInRead, int PosInReadRev, double thisPMD_CT, double thisPMD_GA, BaseContext & Context, int & ReadGroup){throw "Function 'add' Not implemented for base class TSite!"; };
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
	void addToBaseFrequencies(TBaseFrequencies & frequencies);
	double makePhred(double epsilon){
		return makePhredByRef(epsilon);
	};
	double makePhredByRef(double & epsilon){
		if(epsilon < maxQualToPrintNaturalScale) return maxQualToPrint;
		return -10.0 * log10(epsilon);
	};
	void calcEmissionProbabilities();
	void calculateP_g(double* genotypeProbabilities);
	double calculateWeightedSumOfEmissionProbs(double* weights);
	std::string getBases();
	std::string getEmissionProbs();
	double calculateLogLikelihood(double* genotypeProbabilities);
	//MLE Callers
	void calculateNormalizedGenotypeLikelihoods(TRandomGenerator & randomGenerator, double & quality, double & maxGenotypeProb, int & MLGenotype);
	void callMLEGenotype(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool printRef);
	void callMLEGenotypeVCF(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool printRef);
	virtual void callMLEGenotypeKnownAlleles(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){ throw "callMLEGenotypeKnownAlleles not implemented for TSite base class!";};
	virtual void callMLEGenotypeVCFKnownAlleles(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){ throw "callMLEGenotypeVCFKnownAlleles not implemented for TSite base class!";};
	//Bayesian Callers
	void calculateGenotypePosteriorProbabilities(double* pGenotype, TRandomGenerator & randomGenerator, double* postProb, int & MAP);
	void callBayesianGenotype(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool printRef);
	void callBayesianGenotypeVCF(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out);
	virtual void callBayesianGenotypeKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){ throw "callBayesianGenotypeKnownAlleles not implemented for TSite base class!";};
	virtual void callBayesianGenotypeVCFKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){ throw "callBayesianGenotypeVCFKnownAlleles not implemented for TSite base class!";};
	//Allele Presence Callers
	virtual void callAllelePresence(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool printRef){ throw "callAllelePresence not implemented for TSite base class!";};
	virtual void callAllelePresenceKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){ throw "callAllelePresenceKnownAlleles not implemented for TSite base class!";};
	virtual void callAllelePresenceVCF(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out){ throw "callAllelePresenceVCF not implemented for TSite base class!";};
	virtual void callAllelePresenceVCFKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){ throw "callAllelePresenceVCFKnownAlleles not implemented for TSite base class!";};
};

class TSiteDiploid:public TSite{
public:

	TSiteDiploid(){
		hasData = false;
		numGenotypes = 10;
		emissionProbabilities = new double[numGenotypes];
		P_g = new double[numGenotypes];
	};
	~TSiteDiploid(){
		delete[] emissionProbabilities;
		delete[] P_g;
	};
	void add(char & base, char & quality, int PosInRead, int PosInReadRev, double thisPMD_CT, double thisPMD_GA, BaseContext & Context, int & ReadGroup);
	//MLE callers
	void calculatePhredScaledGenotypeLikelihoodsKnownAlleles(TGenotypeMap & genoMap, char & alt, TRandomGenerator & randomGenerator, double* phredEmissionProbs, double & quality, double & maxGenotypeProb, int & MLGenotype);
	void callMLEGenotypeKnownAlleles(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
	void callMLEGenotypeVCFKnownAlleles(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
	//Bayesian Callers
	void calculateGenotypePosteriorProbabilitiesKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, char & alt, TRandomGenerator & randomGenerator, double* postProb, int & MAP);
	void callBayesianGenotypeKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
	void callBayesianGenotypeVCFKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
	//Allele Presence Callers
	void calculatePosteriorOnAllelePresence(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, double* postProbAllele, int & MAP);
	void callAllelePresence(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool printRef);
	void callAllelePresenceVCF(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out);
	void calculatePosteriorOnAllelePresenceKnownAlleles(double* pGenotype, char & alt, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, double* postProbAllele, int & MAP);
	void callAllelePresenceKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
	void callAllelePresenceVCFKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt);
};

class TSiteHaploid:public TSite{
public:

	TSiteHaploid(){
		hasData = false;
		numGenotypes = 4;
		emissionProbabilities = new double[numGenotypes];
		P_g = new double[numGenotypes];
	}
	~TSiteHaploid(){
		delete[] emissionProbabilities;
		delete[] P_g;
	};
	void add(char & base, char & quality, int PosInRead, int PosInReadRev, double thisPMD_CT, double thisPMD_GA, BaseContext & Context, int & ReadGroup);
};


#endif /* TSITE_H_ */
