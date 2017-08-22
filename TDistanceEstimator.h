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
	//Genotype*** phiToGenoMap;
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
//TDistanceEstimate
//--------------------------------------------
class TDistanceEstimate{
private:
	TGenotypeMap genoMap;
	TGenoToPhiMap genoToPhiMap;
	TPhredToLikelihood phredToLik;
	TBaseFrequencies pi1;
	TBaseFrequencies pi2;
	TBaseFrequencies* old_pi;
	TBaseFrequencies* tmp_pi;
	double* old_phi;
	double* tmp_phi;
	double old_LL;
	double** probGeno;
	double* K; //normalizing constant for W
	int g1, g2; //index variables

	void fill_K(TBaseFrequencies  & thesePi);
	void fill_P_g_given_phi_pi(double* phi, TBaseFrequencies & pi);
	void estimatePhiWithEM(int** genoQual1, int** genoQual2, long numSites, int maxNumIterations);

public:
	TBaseFrequencies* pi;
	double* phi;
	double LL;

	TDistanceEstimate();
	~TDistanceEstimate(){
		delete[] phi;
		delete[] old_phi;
		delete[] K;
		for(g1=0; g1<10; ++g1)
			delete[] probGeno[g1];
		delete[] probGeno;
	};
	void guessPi(int** genoQual1, int** genoQual2, long numSites);
	void guessPhi(int** genoQual1, int** genoQual2, long numSites);
};

//--------------------------------------------
//TDistanceEstimator
//--------------------------------------------
class TDistanceEstimator{
private:
	TLog* logfile;

	void estimateDistanceInWindows(std::string filename, TGlfReader & g1, TGlfReader & g2, long windowLen);

	void writeDistanceEstimates(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd, int & numsitesWithData);
	void writeDistanceEstimatesNoData(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd);

public:
	TDistanceEstimator(TLog* Logfile);
	~TDistanceEstimator(){};

	void printGLF(TParameters & params);
	void estimateDistances(TParameters & params);

};


#endif /* TDISTANCEESTIMATOR_H_ */
