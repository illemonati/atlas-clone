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
	void setRefBase(char & Base){referenceBase = Base; };
	char getRefBase(){return referenceBase;};
	Base getRefBaseAsEnum(){
		if(referenceBase == 'A') return A;
		if(referenceBase == 'C') return C;
		if(referenceBase == 'G') return G;
		if(referenceBase == 'T') return T;
		return N;
	};
	void addToBaseFrequencies(TBaseFrequencies & frequencies);
	void calcEmissionProbabilities();
	void callMLEGenotype(TGenotypeMap & genoMap, gz::ogzstream & out, bool printRef);
	virtual void callBayesianGenotype(double* pGenotype, TGenotypeMap & genoMap, gz::ogzstream & out, bool printRef);
	virtual void callAllelePresence(double* pGenotype, TGenotypeMap & genoMap, gz::ogzstream & out, bool printRef){ throw "callMLEAllelePresence not implemented for TSite base class!";};
	void calculateP_g(double* genotypeProbabilities);
	double calculateWeightedSumOfEmissionProbs(double* weights);
	std::string getBases();
	std::string getEmissionProbs();
	double calculateLogLikelihood(double* genotypeProbabilities);
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
	void callAllelePresence(double* pGenotype, TGenotypeMap & genoMap, gz::ogzstream & out, bool printRef);
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
