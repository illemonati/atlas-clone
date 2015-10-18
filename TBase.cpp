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

void TBaseDiploidA::fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbGA(pos3);

	emissionProbabilities.set(AA, oneMinusError);
	emissionProbabilities.set(CC, errorOneThird);
	emissionProbabilities.set(GG, (1.0 - pmd) * errorOneThird + pmd * oneMinusError);
	emissionProbabilities.set(TT, errorOneThird);
	emissionProbabilities.set(AC, 0.5 - errorOneThird);
	emissionProbabilities.set(AG, ((1.0 + pmd) * oneMinusError + (1.0 - pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(AT, 0.5 - errorOneThird);
	emissionProbabilities.set(CG, (pmd * oneMinusError + (2.0 - pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(CT, errorOneThird);
	emissionProbabilities.set(GT, emissionProbabilities.get(CG));
};

void TBaseHaploidA::fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate){
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbGA(pos3);

	emissionProbabilities.set(A, oneMinusError);
	emissionProbabilities.set(C, errorOneThird);
	emissionProbabilities.set(G, pmd * oneMinusError + (1.0 - pmd) * errorOneThird);
	emissionProbabilities.set(T, errorOneThird);
}

void TBaseDiploidC::fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbCT(pos5);

	emissionProbabilities.set(AA , errorOneThird);
	emissionProbabilities.set(CC , (1.0 - pmd) * oneMinusError + pmd * errorOneThird);
	emissionProbabilities.set(GG , errorOneThird);
	emissionProbabilities.set(TT , errorOneThird);
	emissionProbabilities.set(AC , ((1.0 - pmd) * oneMinusError + (1.0 + pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(AG , errorOneThird);
	emissionProbabilities.set(AT , errorOneThird);
	emissionProbabilities.set(CG , emissionProbabilities.get(AC));
	emissionProbabilities.set(CT , emissionProbabilities.get(AC));
	emissionProbabilities.set(GT , errorOneThird);
};

void TBaseHaploidC::fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbCT(pos5);

	emissionProbabilities.set(A, errorOneThird);
	emissionProbabilities.set(C, (1.0 - pmd) * oneMinusError + pmd * errorOneThird);
	emissionProbabilities.set(G, errorOneThird);
	emissionProbabilities.set(T, errorOneThird);
}

void TBaseDiploidG::fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbGA(pos3);

	emissionProbabilities.set(AA, errorOneThird);
	emissionProbabilities.set(CC, errorOneThird);
	emissionProbabilities.set(GG, (1.0 - pmd) * oneMinusError + pmd * errorOneThird);
	emissionProbabilities.set(TT, errorOneThird);
	emissionProbabilities.set(AC, errorOneThird);
	emissionProbabilities.set(AG, ((1.0 - pmd) * oneMinusError + (1.0 + pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(AT, errorOneThird);
	emissionProbabilities.set(CG, emissionProbabilities.get(AG));
	emissionProbabilities.set(CT, errorOneThird);
	emissionProbabilities.set(GT, emissionProbabilities.get(AG));

};

void TBaseHaploidG::fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbGA(pos3);

	emissionProbabilities.set(A, errorOneThird);
	emissionProbabilities.set(C, errorOneThird);
	emissionProbabilities.set(G, (1.0 - pmd) * oneMinusError + pmd * errorOneThird);
	emissionProbabilities.set(T, errorOneThird);
}

void TBaseDiploidT::fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbCT(pos5);

	emissionProbabilities.set(AA, errorOneThird);
	emissionProbabilities.set(CC, (1.0 - pmd) * errorOneThird + pmd * oneMinusError);
	emissionProbabilities.set(GG, errorOneThird);
	emissionProbabilities.set(TT, oneMinusError);
	emissionProbabilities.set(AC, (pmd * oneMinusError + (2.0 - pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(AG, errorOneThird);
	emissionProbabilities.set(AT, 0.5 - errorOneThird);
	emissionProbabilities.set(CG, emissionProbabilities.get(AC));
	emissionProbabilities.set(CT, ((1.0 + pmd) * oneMinusError + (1.0 - pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(GT, 0.5 - errorOneThird);
};

void TBaseHaploidT::fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;
	double pmd = pmdObject.getProbCT(pos5);

	emissionProbabilities.set(A, errorOneThird);
	emissionProbabilities.set(C, pmd * oneMinusError + (1.0 - pmd) * errorOneThird);
	emissionProbabilities.set(G, errorOneThird);
	emissionProbabilities.set(T, oneMinusError);
}


