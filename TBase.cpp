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

void TBase::addToEmissionProbProduct(double* vec){
	double errorOneThird = errorRate / 3.0;
	double oneMinusError = 1.0 - errorRate;

	if(base == A){
		vec[AA] *= oneMinusError;
		vec[CC] *= errorOneThird;
		vec[GG] *= (1.0 - PMD_GA) * errorOneThird + PMD_GA * oneMinusError;
		vec[TT] *= errorOneThird;
		vec[AC] *= 0.5 - errorOneThird;
		vec[AG] *= ((1.0 + PMD_GA) * oneMinusError + (1.0 - PMD_GA) * errorOneThird) / 2.0;
		vec[AT] *= 0.5 - errorOneThird;
		vec[CT] *= errorOneThird;
		double tmp = (PMD_GA * oneMinusError + (2.0 - PMD_GA) * errorOneThird) / 2.0;
		vec[CG] *= tmp;
		vec[GT] *= tmp;
	} else if (base == C){
		vec[AA] *= errorOneThird;
		vec[CC] *= (1.0 - PMD_CT) * oneMinusError + PMD_CT * errorOneThird;
		vec[GG] *= errorOneThird;
		vec[TT] *= errorOneThird;
		vec[AG] *= errorOneThird;
		vec[AT] *= errorOneThird;
		vec[GT] *= errorOneThird;
		double tmp = ((1.0 - PMD_CT) * oneMinusError + (1.0 + PMD_CT) * errorOneThird) / 2.0;
		vec[AC] *= tmp;
		vec[CG] *= tmp;
		vec[CT] *= tmp;
	} else if (base == G){
		vec[AA] *= errorOneThird;
		vec[CC] *= errorOneThird;
		vec[GG] *= (1.0 - PMD_GA) * oneMinusError + PMD_GA * errorOneThird;
		vec[TT] *= errorOneThird;
		vec[AC] *= errorOneThird;
		vec[AT] *= errorOneThird;
		vec[CT] *= errorOneThird;
		double tmp = ((1.0 - PMD_GA) * oneMinusError + (1.0 + PMD_GA) * errorOneThird) / 2.0;
		vec[AG] *= tmp;
		vec[CG] *= tmp;
		vec[GT] *= tmp;
	} else if (base == T){
		vec[AA] *= errorOneThird;
		vec[CC] *= (1.0 - PMD_CT) * errorOneThird + PMD_CT * oneMinusError;
		vec[GG] *= errorOneThird;
		vec[TT] *= oneMinusError;
		vec[AG] *= errorOneThird;
		vec[AT] *= 0.5 - errorOneThird;
		double tmp = (PMD_CT * oneMinusError + (2.0 - PMD_CT) * errorOneThird) / 2.0;
		vec[AC] *= tmp;
		vec[CG] *= tmp;
		vec[CT] *= ((1.0 + PMD_CT) * oneMinusError + (1.0 - PMD_CT) * errorOneThird) / 2.0;
		vec[GT] *= 0.5 - errorOneThird;
	}
}

void TBase::addToEmissionProbSum(double vec[10]){
	double errorOneThird = errorRate / 3.0;
	double oneMinusError = 1.0 - errorRate;

	if(base == A){
		vec[AA] += log(oneMinusError);
		vec[CC] += log(errorOneThird);
		vec[GG] += log((1.0 - PMD_GA) * errorOneThird + PMD_GA * oneMinusError);
		vec[TT] += log(errorOneThird);
		vec[AC] += log(0.5 - errorOneThird);
		vec[AG] += log(((1.0 + PMD_GA) * oneMinusError + (1.0 - PMD_GA) * errorOneThird) / 2.0);
		vec[AT] += log(0.5 - errorOneThird);
		vec[CT] += log(errorOneThird);
		double tmp = log((PMD_GA * oneMinusError + (2.0 - PMD_GA) * errorOneThird) / 2.0);
		vec[CG] += tmp;
		vec[GT] += tmp;
	} else if (base == C){
		vec[AA] += log(errorOneThird);
		vec[CC] += log((1.0 - PMD_CT) * oneMinusError + PMD_CT * errorOneThird);
		vec[GG] += log(errorOneThird);
		vec[TT] += log(errorOneThird);
		vec[AG] += log(errorOneThird);
		vec[AT] += log(errorOneThird);
		vec[GT] += log(errorOneThird);
		double tmp = log(((1.0 - PMD_CT) * oneMinusError + (1.0 + PMD_CT) * errorOneThird) / 2.0);
		vec[AC] += tmp;
		vec[CG] += tmp;
		vec[CT] += tmp;
	} else if (base == G){
		vec[AA] += log(errorOneThird);
		vec[CC] += log(errorOneThird);
		vec[GG] += log((1.0 - PMD_GA) * oneMinusError + PMD_GA * errorOneThird);
		vec[TT] += log(errorOneThird);
		vec[AC] += log(errorOneThird);
		vec[AT] += log(errorOneThird);
		vec[CT] += log(errorOneThird);
		double tmp = ((1.0 - PMD_GA) * oneMinusError + (1.0 + PMD_GA) * errorOneThird) / 2.0;
		vec[AG] += log(tmp);
		vec[CG] += log(tmp);
		vec[GT] += log(tmp);
	} else if (base == T){
		vec[AA] += log(errorOneThird);
		vec[CC] += log((1.0 - PMD_CT) * errorOneThird + PMD_CT * oneMinusError);
		vec[GG] += log(errorOneThird);
		vec[TT] += log(oneMinusError);
		vec[AG] += log(errorOneThird);
		vec[AT] += log(0.5 - errorOneThird);
		double tmp = log(PMD_CT * oneMinusError + (2.0 - PMD_CT) * errorOneThird) / 2.0;
		vec[AC] += tmp;
		vec[CG] += tmp;
		vec[CT] += log(((1.0 + PMD_CT) * oneMinusError + (1.0 - PMD_CT) * errorOneThird) / 2.0);
		vec[GT] += log(0.5 - errorOneThird);
	}
}
/*
void TBase::fillEmissionProbabilitiesCore(double thisErrorRate){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = thisErrorRate / 3.0;
	double oneMinusError = 1.0 - thisErrorRate;

	if(base == A){
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
	} else if (base == C){
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
	} else if (base == G){
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
	} else if (base == T){
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
	}
};
*/


