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


int TSimulatorQuality::returnQual(int qual, int pos, BaseContext baseContext, int maxQual){
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
		tmp = pow(10.0, -(double) (i - 33.0) / 10.0);
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

int TSimulatorRecalTransform::returnQual(int qual, int pos, BaseContext baseContext, int maxQual){
	return transformQuality(qual, pos, baseContext);
}


//--------------------------
//BQSR transformation
//-------------------------

TSimulatorBQSRTransform::TSimulatorBQSRTransform(std::string QualTransform, TSimulatorReadLength* ReadLengthDist): TSimulatorQuality(ReadLengthDist){
	//parse qualTransform input
	std::string::size_type pos = QualTransform.find_first_of('[');
	if(pos == std::string::npos) throw "Can not initialize quality transformation function '" + QualTransform + "': wrong format!\n" + "logit[alpha,beta]";
	std::string name = QualTransform.substr(0,pos);
	if(name == "logit"){
		//prepare first value
		std::string tmp = QualTransform.substr(pos+1, QualTransform.length() - pos - 1);
		pos = tmp.find_first_of(',');
		if(pos == std::string::npos) throw "Can not initialize quality transformation function '" + QualTransform + "': wrong format!\n" + "logit[alpha,beta]";
		alpha = atof(tmp.substr(0, pos).c_str());
		beta = atof(tmp.substr(pos+1).c_str());
	}
}

int TSimulatorBQSRTransform::returnQual(int qual, int pos, BaseContext baseContext, int maxQual){
	double trueError = pow(10,qual/-10);
	double inverseLogit = log((1-trueError)/trueError);
	double fakeError = (1/(1+exp(alpha*inverseLogit*inverseLogit + beta*inverseLogit))); //logit function -> natural logarithm!
	if(fakeError == 0) return maxQual;
	qual = -10 * log10(fakeError);
	return qual;
}

//--------------------------
//BQSR transformation Position
//-------------------------

TSimulatorBQSRPositionTransform::TSimulatorBQSRPositionTransform(float PositionTransform, std::string QualTransform, TSimulatorReadLength* ReadLengthDist): TSimulatorBQSRTransform(QualTransform, ReadLengthDist){
	revIntercept = PositionTransform;
	findTransformationSlope();
}

float TSimulatorBQSRPositionTransform::calculateSum(float & curSlope, float & curIntercept, float & newSum){
	//calculate sum( prob of observing given position) * log(positionBeta) -> should be 0
	newSum = 1.0;
	for(int p=0; p<readLengthDist->max(); ++p){
//		std::cout << "curIntercept: " << curIntercept << std::endl;
//		std::cout << log(curSlope * p + curIntercept) << std::endl;
		newSum *= (1.0 - readLengthDist->gammaCumulDensity[p]) * log(curSlope * p + curIntercept);
//		std::cout << "newSum = " << newSum << ", curSlope: " << curSlope << "; curIntercept: " << curIntercept << std::endl ;
		if(curIntercept < 0) throw "curIntercept is smaller than 0!";
		if((curSlope * p + curIntercept) < 0) throw "beta position is negative! curSlope * p + curIntercept: " + toString(curSlope * p + curIntercept) + ", curSlope: " + toString(curSlope) + " p: " + toString(p) + " curIntercept: " + toString(curIntercept);
	}
	return newSum;
}

void TSimulatorBQSRPositionTransform::findTransformationSlope(){
	int numIterations = 25;
	int x = readLengthDist->max();
	float y = revIntercept;
	float curSlope = 0.001, curStepSize = 0.2, newSum = -1.0;
	float curIntercept = -curSlope * x + y;

	float curSum = calculateSum(curSlope, curIntercept, newSum);

	for(int i=0; i<numIterations; ++i){
		std::cout << "######### i " << i << std::endl;
		if(curSum == 0){
			std::cout << "curSum = 0" << std::endl;
			break;
		}

		std::cout << "stepSize: " << curStepSize << " curIntercept: " << curIntercept << " curSum: " << curSum << " curSlope: " << curSlope << std::endl;
		float a = -(curSlope + curStepSize) * x + y;
		std::cout << "(-curSlope + curStepSize) * x + y): " << a << std::endl;
		while((curSlope + curStepSize) < 0.0 || a < 0.0){ //check that neither slope nor intercept become negative
			std::cout << "in while loop" << std::endl;
			curStepSize = curStepSize / std::exp(1.0);
			if(curStepSize < 0.00000000000000001) break;
		}
		curSlope = curSlope + curStepSize;
		curSum = newSum;
		curIntercept = -curSlope * x + y;

		std::cout << "stepSize: " << curStepSize << " curIntercept: " << curIntercept << " curSum: " << curSum << " curSlope: " << curSlope << std::endl;

		newSum = calculateSum(curSlope, curIntercept, newSum);

		if((newSum > 0.0 && curSum < 0.0) || (newSum < 0.0 && curSum > 0.0)){
			std::cout << "changing step direction" << std::endl;
			curStepSize = -curStepSize / std::exp(1.0);
		}
		std::cout << "stepSize: " << curStepSize << " curIntercept: " << curIntercept << " curSum: " << curSum << " newSum: " << newSum << " curSlope: " << curSlope << std::endl;
	}
}




















