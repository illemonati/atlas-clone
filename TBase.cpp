/*
 * TBase.cpp
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */


#include "TBase.h"

//-------------------------------------------------------
//TBase
//-------------------------------------------------------
/*
void TBase::fillEmissionProbabilities(TPMD & pmdObject){
	fillEmissionProbabilitiesCore(pmdObject, errorRate);
}
*/

void TBaseDiploidA::fillEmissionProbabilitiesCore(double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;

	emissionProbabilities.set(AA, oneMinusError);
	emissionProbabilities.set(CC, errorOneThird);
	emissionProbabilities.set(GG, (1.0 - PMD_GA) * errorOneThird + PMD_GA * oneMinusError);
	emissionProbabilities.set(TT, errorOneThird);
	emissionProbabilities.set(AC, 0.5 - errorOneThird);
	emissionProbabilities.set(AG, ((1.0 + PMD_GA) * oneMinusError + (1.0 - PMD_GA) * errorOneThird) / 2.0);
	emissionProbabilities.set(AT, 0.5 - errorOneThird);
	emissionProbabilities.set(CG, (PMD_GA * oneMinusError + (2.0 - PMD_GA) * errorOneThird) / 2.0);
	emissionProbabilities.set(CT, errorOneThird);
	emissionProbabilities.set(GT, emissionProbabilities.get(CG));
};

void TBaseHaploidA::fillEmissionProbabilitiesCore(double thisErrorRate){
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;

	emissionProbabilities.set(A, oneMinusError);
	emissionProbabilities.set(C, errorOneThird);
	emissionProbabilities.set(G, PMD_GA * oneMinusError + (1.0 - PMD_GA) * errorOneThird);
	emissionProbabilities.set(T, errorOneThird);
}

void TBaseDiploidC::fillEmissionProbabilitiesCore(double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;

	emissionProbabilities.set(AA , errorOneThird);
	emissionProbabilities.set(CC , (1.0 - PMD_CT) * oneMinusError + PMD_CT * errorOneThird);
	emissionProbabilities.set(GG , errorOneThird);
	emissionProbabilities.set(TT , errorOneThird);
	emissionProbabilities.set(AC , ((1.0 - PMD_CT) * oneMinusError + (1.0 + PMD_CT) * errorOneThird) / 2.0);
	emissionProbabilities.set(AG , errorOneThird);
	emissionProbabilities.set(AT , errorOneThird);
	emissionProbabilities.set(CG , emissionProbabilities.get(AC));
	emissionProbabilities.set(CT , emissionProbabilities.get(AC));
	emissionProbabilities.set(GT , errorOneThird);
};

void TBaseHaploidC::fillEmissionProbabilitiesCore(double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;

	emissionProbabilities.set(A, errorOneThird);
	emissionProbabilities.set(C, (1.0 - PMD_CT) * oneMinusError + PMD_CT * errorOneThird);
	emissionProbabilities.set(G, errorOneThird);
	emissionProbabilities.set(T, errorOneThird);
}

void TBaseDiploidG::fillEmissionProbabilitiesCore(double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;

	emissionProbabilities.set(AA, errorOneThird);
	emissionProbabilities.set(CC, errorOneThird);
	emissionProbabilities.set(GG, (1.0 - PMD_GA) * oneMinusError + PMD_GA * errorOneThird);
	emissionProbabilities.set(TT, errorOneThird);
	emissionProbabilities.set(AC, errorOneThird);
	emissionProbabilities.set(AG, ((1.0 - PMD_GA) * oneMinusError + (1.0 + PMD_GA) * errorOneThird) / 2.0);
	emissionProbabilities.set(AT, errorOneThird);
	emissionProbabilities.set(CG, emissionProbabilities.get(AG));
	emissionProbabilities.set(CT, errorOneThird);
	emissionProbabilities.set(GT, emissionProbabilities.get(AG));

};

void TBaseHaploidG::fillEmissionProbabilitiesCore(double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;

	emissionProbabilities.set(A, errorOneThird);
	emissionProbabilities.set(C, errorOneThird);
	emissionProbabilities.set(G, (1.0 - PMD_GA) * oneMinusError + PMD_GA * errorOneThird);
	emissionProbabilities.set(T, errorOneThird);
}

void TBaseDiploidT::fillEmissionProbabilitiesCore(double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;

	emissionProbabilities.set(AA, errorOneThird);
	emissionProbabilities.set(CC, (1.0 - PMD_CT) * errorOneThird + PMD_CT * oneMinusError);
	emissionProbabilities.set(GG, errorOneThird);
	emissionProbabilities.set(TT, oneMinusError);
	emissionProbabilities.set(AC, (PMD_CT * oneMinusError + (2.0 - PMD_CT) * errorOneThird) / 2.0);
	emissionProbabilities.set(AG, errorOneThird);
	emissionProbabilities.set(AT, 0.5 - errorOneThird);
	emissionProbabilities.set(CG, emissionProbabilities.get(AC));
	emissionProbabilities.set(CT, ((1.0 + PMD_CT) * oneMinusError + (1.0 - PMD_CT) * errorOneThird) / 2.0);
	emissionProbabilities.set(GT, 0.5 - errorOneThird);
};

void TBaseHaploidT::fillEmissionProbabilitiesCore(double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;

	emissionProbabilities.set(A, errorOneThird);
	emissionProbabilities.set(C, PMD_CT * oneMinusError + (1.0 - PMD_CT) * errorOneThird);
	emissionProbabilities.set(G, errorOneThird);
	emissionProbabilities.set(T, oneMinusError);
}


