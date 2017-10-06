/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#include "TSimulatorQuality.h"

TSimulatorQuality::TSimulatorQuality(TSimulatorReadLength* ReadLengthDist){
	readLengthDist = ReadLengthDist;
};


int TSimulatorQuality::returnQual(int qual, int pos, BaseContext baseContext){
	return qual;
}


//--------------------------
//recal transformation
//-------------------------
TSimulatorRecalTransform::TSimulatorRecalTransform(std::vector<double> Betas, TSimulatorReadLength* ReadLengthDist): TSimulatorQuality(ReadLengthDist){
	if(Betas.size() != 24)
		throw "Wrong size of beta vector when initializing quality transformation: need 24 values (quality, quality^2, pos, pos^2 and 20 contexts).";

	//copy betas
	beta = new double[24];
	for(int i=0; i<24; ++i)
		beta[i] = Betas[i];

	//precalculate stuff
	qualTermForTransformation = new double[127];
	for(int i=0; i<34; ++i)
		qualTermForTransformation[i] = 100.0;

	double tmp;
	for(int i=34; i<127; ++i){
		tmp = pow(10.0, -(double) (i - 33) / 10.0);
		qualTermForTransformation[i] = log(tmp / (1.0 - tmp));
	}
	posTermForTransformation = new double[readLengthDist->max()];

	for(int i=0; i<readLengthDist->max(); ++i){
		posTermForTransformation[i] = beta[2] * i + beta[3] * i*i;
	}
}

int TSimulatorRecalTransform::transformQuality(int & qual, int pos, int context){
	static double constant;
	static double tmp;
	static double q;

	constant = posTermForTransformation[pos] + beta[context+4] - qualTermForTransformation[qual];
	if(4.0 * beta[1] * constant > beta[0] * beta[0]) throw "beta[0]^2 cannot be smaller than 4beta[1](position + context constants)";
	if(beta[1] == 0.0){
		q = -constant / beta[0];
	} else {
		tmp = sqrt(beta[0] * beta[0] - 4.0 * beta[1] * constant);
		q = (-tmp - beta[0]) / 2.0 / beta[1];
		//if(q < 0) q = (-tmp - beta[0]) / 2.0 / beta[1];
	}

	tmp = exp(q);
	if(tmp == 0) throw "choose different quality transformation parameters! tmp == 0";
	return -10.0 * log10(tmp / (1.0 + tmp)) + 33.0;
}

int TSimulatorRecalTransform::returnQual(int qual, int pos, BaseContext baseContext){
	return transformQuality(qual, pos, baseContext);
}
