/*
 * TDistanceCalculator.h
 *
 *  Created on: Jul 25, 2017
 *      Author: phaentu
 */

#ifndef TDISTANCEESTIMATOR_H_
#define TDISTANCEESTIMATOR_H_

#include "TParameters.h"
#include "TGLF.h"
#include <math.h>
#include "TGenotypeMap.h"

//--------------------------------------------
//TPhiToGenoMap
//--------------------------------------------
class TGenoToPhiMap{
public:
	TGenotypeMap genoMap;
	int**   genoToPhiMap;

	TGenoToPhiMap();
	~TGenoToPhiMap(){
		for(int i=0; i<10; ++i)
			delete[] genoToPhiMap[i];
		delete[] genoToPhiMap;
	};

	int& operator()(int & g1, int & g2){
		return genoToPhiMap[g1][g2];
	};
	int& operator()(Genotype & g1, Genotype & g2){
		return genoToPhiMap[g1][g2];
	};
};

//--------------------------------------------
//TGenocombinationToBaseMap
//--------------------------------------------
class TGenocombinationToBaseMap{
public:
	TGenotypeMap genoMap;
	bool*** genotypeCombinationHasBase;

	TGenocombinationToBaseMap();
	~TGenocombinationToBaseMap(){
		for(int i=0; i<10; ++i){
			for(int j=0; j<10; ++j)
				delete[] genotypeCombinationHasBase[i][j];
		}
		delete[] genotypeCombinationHasBase;
	};

	bool& operator()(int & g1, int & g2, int & base){
		return genotypeCombinationHasBase[g1][g2][base];
	};
	bool& operator()(Genotype & g1, Genotype & g2, Base & b){
		return genotypeCombinationHasBase[g1][g2][b];
	};
};

//--------------------------------------------
//TDistanceEstimate
//--------------------------------------------
class TEMforDistanceEstimation{
private:
	TGenotypeMap genoMap;
	TGenoToPhiMap genoToPhiMap;
	TGenocombinationToBaseMap genoToBaseMap;
	TPhredToLikelihood phredToLik;
	double old_LL;
	double* K; //normalizing constant
	double** probGeno;
	double** P_G;
	double** P_G_one_site;
	int g1, g2; //index variables

	void guessPi(int** genoQual1, int** genoQual2, long numSites);
	void guessPhi(int** genoQual1, int** genoQual2, long numSites);
	void fill_K(TBaseFrequencies  & thesePi);
	void fill_P_g_given_phi_pi(double* phi, TBaseFrequencies & pi);

public:
	TBaseFrequencies pi;
	double* phi;
	double LL;

	TEMforDistanceEstimation();
	~TEMforDistanceEstimation(){
		delete[] phi;
		delete[] K;
		for(g1=0; g1<10; ++g1){
			delete[] probGeno[g1];
			delete[] P_G[g1];
			delete[] P_G_one_site[g1];
		}
		delete[] probGeno;
		delete[] P_G;
		delete[] P_G_one_site;
	};

	bool estimatePhiWithEM(int** genoQual1, int** genoQual2, long numSites, int maxNumIterations, double epsilon);
};

//--------------------------------------------
//TDistanceEstimator
//--------------------------------------------
class TDistanceEstimator{
private:
	TLog* logfile;
	int maxNumEMIterations;
	double epsilonForEM;

	void estimateDistanceInWindows(std::string filename, TGlfReader & g1, TGlfReader & g2, long windowLen);

	void writeDistanceEstimates(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd, int & numsitesWithData, TEMforDistanceEstimation & EM_object);
	void writeDistanceEstimatesNoData(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd);

public:
	TDistanceEstimator(TLog* Logfile);
	~TDistanceEstimator(){};

	void printGLF(TParameters & params);
	void estimateDistances(TParameters & params);

};


#endif /* TDISTANCEESTIMATOR_H_ */
