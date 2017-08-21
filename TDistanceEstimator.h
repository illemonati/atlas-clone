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
class TPhiToGenoMap{
public:
	//Genotype*** phiToGenoMap;
	int** genoToPhiMap;

	TPhiToGenoMap(){
		genoToPhiMap = new int*[10];
		for(int i=0; i<10; ++i)
			genoToPhiMap[i] = new int[10];

		//case aa/aa
		genoToPhiMap[AA][AA] = 0;
		genoToPhiMap[CC][CC] = 0;
		genoToPhiMap[GG][GG] = 0;
		genoToPhiMap[TT][TT] = 0;

		//case aa/ab
		genoToPhiMap[AA][AC] = 1;
		genoToPhiMap[AA][AG] = 1;
		genoToPhiMap[AA][AT] = 1;
		genoToPhiMap[CC][AC] = 1;
		genoToPhiMap[CC][CG] = 1;
		genoToPhiMap[CC][CT] = 1;
		genoToPhiMap[GG][AG] = 1;
		genoToPhiMap[GG][CG] = 1;
		genoToPhiMap[GG][GT] = 1;
		genoToPhiMap[TT][AT] = 1;
		genoToPhiMap[TT][CT] = 1;
		genoToPhiMap[TT][GT] = 1;

		//case ab/aa
		genoToPhiMap[AC][AA] = 2;
		genoToPhiMap[AG][AA] = 2;
		genoToPhiMap[AT][AA] = 2;
		genoToPhiMap[AC][CC] = 2;
		genoToPhiMap[CG][CC] = 2;
		genoToPhiMap[CT][CC] = 2;
		genoToPhiMap[AG][GG] = 2;
		genoToPhiMap[CG][GG] = 2;
		genoToPhiMap[GT][GG] = 2;
		genoToPhiMap[AT][TT] = 2;
		genoToPhiMap[CT][TT] = 2;
		genoToPhiMap[GT][TT] = 2;

		//case aa/bb
		genoToPhiMap[AA][CC] = 3;
		genoToPhiMap[AA][GG] = 3;
		genoToPhiMap[AA][TT] = 3;
		genoToPhiMap[CC][AA] = 3;
		genoToPhiMap[CC][GG] = 3;
		genoToPhiMap[CC][TT] = 3;
		genoToPhiMap[GG][AA] = 3;
		genoToPhiMap[GG][CC] = 3;
		genoToPhiMap[GG][TT] = 3;
		genoToPhiMap[TT][AA] = 3;
		genoToPhiMap[TT][CC] = 3;
		genoToPhiMap[TT][GG] = 3;

		//case ab/ab
		genoToPhiMap[AC][AC] = 4;
		genoToPhiMap[AG][AG] = 4;
		genoToPhiMap[AT][AT] = 4;
		genoToPhiMap[CG][CG] = 4;
		genoToPhiMap[CT][CT] = 4;
		genoToPhiMap[GT][GT] = 4;

		//case ab/cc
		genoToPhiMap[AC][GG] = 5;
		genoToPhiMap[AC][TT] = 5;
		genoToPhiMap[AG][CC] = 5;
		genoToPhiMap[AG][TT] = 5;
		genoToPhiMap[AT][CC] = 5;
		genoToPhiMap[AT][GG] = 5;
		genoToPhiMap[CG][AA] = 5;
		genoToPhiMap[CG][TT] = 5;
		genoToPhiMap[CT][AA] = 5;
		genoToPhiMap[CT][GG] = 5;
		genoToPhiMap[GT][AA] = 5;
		genoToPhiMap[GT][CC] = 5;

		//case aa/bc
		genoToPhiMap[GG][AC] = 6;
		genoToPhiMap[TT][AC] = 6;
		genoToPhiMap[CC][AG] = 6;
		genoToPhiMap[TT][AG] = 6;
		genoToPhiMap[CC][AT] = 6;
		genoToPhiMap[GG][AT] = 6;
		genoToPhiMap[AA][CG] = 6;
		genoToPhiMap[TT][CG] = 6;
		genoToPhiMap[AA][CT] = 6;
		genoToPhiMap[GG][CT] = 6;
		genoToPhiMap[AA][GT] = 6;
		genoToPhiMap[CC][GT] = 6;

		//case ab/ac
		genoToPhiMap[AC][AG] = 7;
		genoToPhiMap[AC][AT] = 7;
		genoToPhiMap[AG][AC] = 7;
		genoToPhiMap[AG][AT] = 7;
		genoToPhiMap[AT][AC] = 7;
		genoToPhiMap[AT][AG] = 7;
		genoToPhiMap[CG][AC] = 7;
		genoToPhiMap[CG][CT] = 7;


		//case ab/cd
		genoToPhiMap[AC][GT] = 8;
		genoToPhiMap[AG][CT] = 8;
		genoToPhiMap[AT][CG] = 8;
		genoToPhiMap[CG][AT] = 8;
		genoToPhiMap[CT][AG] = 8;
		genoToPhiMap[GT][AC] = 8;

		/*
		phiToGenoMap = new Genotype**[9];

		//Case aa/aa
		phiToGenoMap[0] = new Genotype*[];
		phiToGenoMap[0][0] = new Genotype[2];
		phiToGenoMap[0][0][0] = AA;
		*/

	};
	~TPhiToGenoMap(){

	};
};


//--------------------------------------------
//TDistanceEstimate
//--------------------------------------------
class TDistanceEstimate{
private:
	TGenotypeMap genoMap;
	TPhredToLikelihood phredToLik;

public:
	TBaseFrequencies pi, old_pi;
	double* phi;
	double* old_phi;
	double LL;
	double old_LL;

	TDistanceEstimate();
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
