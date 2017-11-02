/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#include "TSimulatorRead.h"
#include "TParameters.h"

std::string TSimulatorRead::getQueryBases(){
	return(queryBases);
}
std::string TSimulatorRead::getBamQualities(){
	return(bamQualities);
}

double TSimulatorRead::dePhred(double x){
	return pow(10, -(x-33.0) / 10.0);
}


/*void TSimulatorRead::setPMD(TPMD* PmdObject){
	pmdObject = PmdObject;
	pmdInitialized = true;
}*/

void TSimulatorRead::setQualityDistribution(double mean, double sd, int maxQ){
	meanQual = mean + 33.0; //add 33 to mean quality to get in in char
	sdQual = sd;
	maxQual = maxQ;
}

int TSimulatorRead::sampleQuality(){
	int qual = round(randomGenerator->getNormalRandom(meanQual, sdQual));
	if(qual > 126) qual = 126;
	if(qual < 33) qual = 33;
	return qual;
};

void TSimulatorRead::initializeQualToErrorTable(){
	if(!qualToErroTableInitialized){
		qualToErroTable = new double[127];
		for(int i=0; i<33; ++i)
			qualToErroTable[i] = 1.0;
		for(int i=33; i<127; ++i)
			qualToErroTable[i] = dePhred(i);
	}
	qualToErroTableInitialized = true;
};

TSimulatorRead::TSimulatorRead(TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator){
	logfile = Logfile;
	randomGenerator = RandomGenerator;

	double mQ = params.getParameterDoubleWithDefault("meanQual", 30);
	double sdQ = params.getParameterDoubleWithDefault("sdQual", 10);
	int maxQ = params.getParameterDoubleWithDefault("maxQual", 500);

	setQualityDistribution(mQ, sdQ, maxQ);
	logfile->list("Will simulate normal distributed quality scores with mean = " + toString(meanQual) + " and sd = " + toString(sdQual) + ", capped at " + toString(maxQual) + ".");

	//PMD
	if(params.parameterExists("pmd") || params.parameterExists("pmdCT") || params.parameterExists("pmdGA")){
		pmdInitialized = true;
		pmdObject = new TPMD(params, logfile);
	} else pmdInitialized = false;

	initializeQualToErrorTable();

};


int TSimulatorRead::returnQual(int qual, int pos, BaseContext baseContext, int maxQual){
	return qual;
}

void TSimulatorRead::simulate(short* posAddress, readLengthContainer & rl){
	static short base;
	static int qual;
	static long p;

	//simulate a read starting here
	std::string queryBases = "";
	std::string qualities = "";
	//choose haplotype
	static int previousBase = 4; //means N
	for(p=0; p<rl.readLength; ++p){
		//get true nucleotide
		base = *(posAddress + p);

		//apply PMD
		if(pmdInitialized){
			if(base == 1 ){ //means is C
				if(randomGenerator->getRand() < pmdObject->getProbCT(p))
					base = 3; //means T
			} else if(base == 2){ //means is G
				if(randomGenerator->getRand() < pmdObject->getProbGA(rl.fragmentLength - p - 1))
					base = 0; //means A
			}
		}
		//sample quality and add error
		qual = sampleQuality();
		if(randomGenerator->getRand() < qualToErroTable[qual])
			base = (base + randomGenerator->pickOne(3) + 1) % 4;

		//add to bam alignment
		//int returnQual(int & qual, int & pos, TGenotypeMap & genoMap, int & previousBase, int & base);
		int transQual = returnQual(qual, p, genoMap.getContext(previousBase, base), maxQual);
		if(transQual > maxQual)	qualities += (char) maxQual;
		else qualities += (char) transQual;
		previousBase = base;
		queryBases += toBase[base];
	}

}

//--------------------------
//recal transformation
//-------------------------

TSimulatorReadRecal::TSimulatorReadRecal(std::vector<double> Betas, int & maxReadLen, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator): TSimulatorRead(params, Logfile, RandomGenerator){
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
	posTermForTransformation = new double[maxReadLen];

	for(int i=0; i<maxReadLen; ++i){
		posTermForTransformation[i] = beta[2] * i + beta[3] * i*i;
	}
}

int TSimulatorReadRecal::transformQuality(int & qual, int & pos, int & context){
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

int TSimulatorReadRecal::returnQual(int & qual, int & pos, int & baseContext){
	return transformQuality(qual, pos, baseContext);
}

/*
//--------------------------
//BQSR transformation
//-------------------------
TSimulatorBQSRTransform::TSimulatorBQSRTransform(int & maxReadLen, TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator): TSimulatorRead(Params, Logfile, RandomGenerator){
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



void TSimulatorBQSRPositionTransform::findTransformationSlope(){
}*/




















