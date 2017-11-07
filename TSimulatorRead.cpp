/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#include "TSimulatorRead.h"
#include "TParameters.h"

/*double TSimulatorRead::dePhredAscii(int x){
	return pow(10, -((double) x-33.0) / 10.0);
}*/

double TSimulatorRead::dePhred(int x){
	return pow(10, -((double) x) / 10.0);
}

int TSimulatorRead::phred(double x){
	return(round(-10.0 * log10(x)));
}

void TSimulatorRead::initializeDePhredTable(){
	if(!dePhredTableInitialized){
		dePhredTable = new double[maxQualPlusOne]; //-minQual
		for(int i=0; i<maxQualPlusOne; ++i)
			dePhredTable[i] = dePhred(i);
	}
	dePhredTableInitialized = true;
};

TSimulatorRead::TSimulatorRead(TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase){
	logfile = Logfile;
	randomGenerator = RandomGenerator;
	toBase = ToBase;

	//PMD
	if(params.parameterExists("pmd") || params.parameterExists("pmdCT") || params.parameterExists("pmdGA")){
		pmdInitialized = true;
		pmdObject = new TPMD(params, logfile);
	} else pmdInitialized = false;

	//qual params
	double mQ = params.getParameterDoubleWithDefault("meanQual", 30.0);
	double sdQ = params.getParameterDoubleWithDefault("sdQual", 10.0);
	int maxQ = params.getParameterIntWithDefault("maxQual", 93);

	setTrueQualityDistribution(mQ, sdQ, maxQ);
	logfile->list("Will simulate normal distributed quality scores with mean = " + toString(meanQual) + " and sd = " + toString(sdQual) + ", capped at " + toString(maxQual) + ".");

	initializeDePhredTable();
	initializeQualToErrorTable();
};

std::string TSimulatorRead::getQueryBases(){
	return(queryBases);
}
std::string TSimulatorRead::getBamQualities(){
	return(bamQualities);
}

void TSimulatorRead::setTrueQualityDistribution(double mean, double sd, int maxQ){
	meanQual = mean; //add 33 later to get quality to get in in char
	sdQual = sd;
	maxQual = maxQ;
	maxQualPlusOne = maxQual + 1;
}

int TSimulatorRead::sampleTrueQuality(){
	int qual = round(randomGenerator->getNormalRandom(meanQual, sdQual));
	if(qual > 93) qual = 93;
	if(qual < 0) qual = 0;
	return qual;
};

void TSimulatorRead::initializeQualToErrorTable(){
	if(!dePhredTableInitialized) throw("Cannot initialize qualToErrorTable without dePhredTable!");
	if(!qualToErroTableInitialized){
		qualToErroTable = new double[maxQualPlusOne];
		for(int i=0; i<maxQualPlusOne; ++i){
			qualToErroTable[i] = 1; //dePhredTable[i];
		}
	}
	qualToErroTableInitialized = true;
};

int TSimulatorRead::returnBamQual(int qual, int pos, BaseContext baseContext, int maxQual){
	return qual;
}

void TSimulatorRead::applyPMD(short & base, long & posInRead, readLengthContainer & rl){
	if(base == 1 ){ //means is C
		if(randomGenerator->getRand() < pmdObject->getProbCT(posInRead)){
			base = 3; //means T
		} else if(base == 2){ //means is G
			if(randomGenerator->getRand() < pmdObject->getProbGA(rl.fragmentLength - posInRead - 1)){
				base = 0; //means A
			}
		}
	}
}

void TSimulatorRead::simulate(short* posAddress, readLengthContainer & rl, TGenotypeMap & genoMap){
	static short base;
	static int qual;
	static long p;
	queryBases = "";
	bamQualities = "";
	static int previousBase = 4; //means N

	for(p=0; p<rl.readLength; ++p){
		//get true nucleotide
		base = *(posAddress + p);

		//apply PMD
		if(pmdInitialized){
			applyPMD(base, p, rl);
		}
		//sample quality and add error
		qual = sampleTrueQuality();
		if(randomGenerator->getRand() < qualToErroTable[qual]){
			base = (base + randomGenerator->pickOne(3) + 1) % 4;
		}
		//add to bam alignment
		int transQual = returnBamQual(qual, p, genoMap.getContext(previousBase, base), maxQual);
		if(transQual > maxQual)	bamQualities += (char) maxQual;
		else bamQualities += (char) transQual + 33;
		previousBase = base;
		queryBases += toBase[base];
	}
}

//--------------------------
//recal transformation
//-------------------------

TSimulatorReadRecal::TSimulatorReadRecal(std::vector<double> Betas, int & maxReadLen, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase): TSimulatorRead(params, Logfile, RandomGenerator, ToBase){
	betas = Betas;
	if(betas.size() != 24)
		throw "Wrong size of beta vector when initializing quality transformation: need 24 values (quality, quality^2, pos, pos^2 and 20 contexts).";

	//precalculate stuff
	qualTermForTransformation = new double[127];
	double tmp;
	for(int i=0; i<93.0; ++i){
		tmp = pow(10.0, -(double) i / 10.0);
		qualTermForTransformation[i] = log(tmp / (1.0 - tmp));
	}
	posTermForTransformation = new double[maxReadLen];

	for(int i=0; i<maxReadLen; ++i){
		posTermForTransformation[i] = betas[2] * i + betas[3] * i*i;
	}
}

int TSimulatorReadRecal::transformQuality(int & qual, int & pos, int & context){
	static double constant;
	static double tmp;
	static double q;

	constant = posTermForTransformation[pos] + betas[context+4] - qualTermForTransformation[qual];
	if(4.0 * betas[1] * constant > betas[0] * betas[0]) throw "beta[0]^2 cannot be smaller than 4beta[1](position + context constants)";
	if(betas[1] == 0.0){
		q = -constant / betas[0];
	} else {
		tmp = sqrt(betas[0] * betas[0] - 4.0 * betas[1] * constant);
		q = (-tmp - betas[0]) / 2.0 / betas[1];
		//if(q < 0) q = (-tmp - beta[0]) / 2.0 / beta[1];
	}

	tmp = exp(q);
	if(tmp == 0) throw "choose different quality transformation parameters! tmp == 0";
	return -10.0 * log10(tmp / (1.0 + tmp));
}

/*int TSimulatorReadRecal::returnQual(int & qual, int & pos, int & baseContext){
	return transformQuality(qual, pos, baseContext);
}*/


//--------------------------
//BQSR base transformation
//-------------------------

TSimulatorReadBQSR::TSimulatorReadBQSR(TSimulatorReadLength* ReadLengthDist, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase): TSimulatorRead(params, Logfile, RandomGenerator, ToBase){
	//position parameters
	readLengthDist = ReadLengthDist;
	maxPos = readLengthDist->max();

	//quality parameters
	parseBQSRQualInput(params);
	setFakeQualityDistribution(); //first find kappa and lambda
	initializeDePhredTable();
	initializeFakeQualToTrueQualTable();
}

void TSimulatorReadBQSR::parseBQSRQualInput(TParameters & params){
	//parse qualTransform input
	std::string transParams = params.getParameterString("BQSRQuality");
	std::string::size_type pos = transParams.find_first_of('[');
	if(pos == std::string::npos) throw "Can not initialize quality transformation function '" + transParams + "': wrong format!\n" + "f[alpha,beta]";
	std::string name = transParams.substr(0,pos);
	if(name == "f"){
		//prepare first value
		std::string tmp = transParams.substr(pos+1, transParams.length() - pos - 1);
		pos = tmp.find_first_of(',');
		if(pos == std::string::npos) throw "Can not initialize quality transformation function '" + transParams + "': wrong format!\n" + "f[alpha,beta]";
		phi1 = atof(tmp.substr(0, pos).c_str());
		phi2 = atof(tmp.substr(pos+1).c_str());
	}

	logfile->list("Simulating a quality distortion with alpha1 = " + toString(phi1) + " and alpha2 = " + toString(phi2));
}

int TSimulatorReadBQSR::sampleFakeQuality(){
	int qual = round(randomGenerator->getNormalRandom(kappa, lambda));
	if(qual > 93) qual = 93;
	if(qual < 0) qual = 0;
	return qual;
}

//---------------------------------
// optimization functions
//---------------------------------

void TSimulatorReadBQSR::fillQBetaQBetaP(){
	int num_of_row = maxPos;
	int num_of_col = maxQualPlusOne - minQual;
	double init_value = -1.0;

	//now we have an empty 2D-matrix of size (0,0). Resizing it with one single command:
	QBetaQBetaP.resize( num_of_col , std::vector<double>( num_of_row , init_value ) );
	for(int q = minQual; q < maxQualPlusOne; ++ q){
		for(int p = 0; p<readLengthDist->max(); ++p){
			QBetaQBetaP[q][p] = phred(dePhredTable[q] * returnBetaPp(p));
		}
	}
	BetaQBetaPInitialized = true;
}

void TSimulatorReadBQSR::fillWeights(double & kappa_cur, double & lambda_cur){
	w = new double[maxQualPlusOne - minQual];

	//w at minQual
	w[0] = randomGenerator->normalCumulativeDistributionFunction(((double) minQual + 0.5), kappa_cur, lambda_cur);

	//w at intermediate Q
	for(int q = (minQual + 1) + 0.5; q < maxQual - 1; ++q){
		double start = randomGenerator->normalCumulativeDistributionFunction((double) q - 0.5, kappa_cur, lambda_cur);
		double end = randomGenerator->normalCumulativeDistributionFunction((double) q + 0.5, kappa_cur, lambda_cur);
		w[q-minQual] = end - start;
	}

	//w at maxQual
	w[(maxQual - minQual)] = 1 - randomGenerator->normalCumulativeDistributionFunction(((double) maxQual - 0.5), kappa_cur, lambda_cur);
	weightsInitialized = true;
}

double TSimulatorReadBQSR::returnCurKappa(){
	float kappa = 0.0;
	for(int q = minQual; q < maxQualPlusOne; ++ q){
		for(int p = 0; p<readLengthDist->max(); ++p){
			kappa += QBetaQBetaP[q][p]* readLengthDist->gammaCumulDensity[p] ; //* w[q]);
		}
	}
	return(kappa);
}

double TSimulatorReadBQSR::returnCurLambda(double & kappa){
	float lambda = 0.0;
	for(int q = minQual; q < maxQualPlusOne; ++ q){
		for(int p = 0; p<readLengthDist->max(); ++p){
			lambda += ((QBetaQBetaP[q][p] - kappa) * (QBetaQBetaP[q][p] - kappa) * readLengthDist->gammaCumulDensity[p] * w[q]);
		}
	}
	return lambda;
}

double TSimulatorReadBQSR::returnDelta(double & kappa, double & lambda){
//	kappa = returnCurKappa();
//	lambda = returnCurLambda();

	delta = (kappa - meanQual)*(kappa - meanQual) + (lambda - sdQual)*(lambda - sdQual);
	return(delta);
}

double TSimulatorReadBQSR::updateParam(double & param, float & stepSize, int & nIter){
	for(int i=0; i<nIter; ++i){
		delta_old = delta;
		param += stepSize;
//		delta = returnDelta();
		if(delta > delta_old)
		delta = -delta/exp(1);
	}
	return(param);
}

//---------------------------------

void TSimulatorReadBQSR::setFakeQualityDistribution(){
	double kappa_cur = meanQual;
	double lambda_cur = sdQual;

	fillQBetaQBetaP();
	fillWeights(kappa_cur, lambda_cur);

	delta = returnDelta(kappa_cur, lambda_cur);

	int nTurns = 10;
	int nIter = 50;

	for(int t=0; t<nTurns; ++t){
		//update kappa
		float stepSize = 5.0;
		kappa_cur = updateParam(kappa_cur, stepSize, nIter);

		//update lambda
		stepSize = 1.0;
		lambda_cur = updateParam(lambda_cur, stepSize, nIter);
	}

	kappa = kappa_cur;
	lambda = lambda_cur;
}

int TSimulatorReadBQSR::returnTrueQual(int & fakeQual){
	double fakeError = dePhredTable[fakeQual];
	return(pow(10, -1/10 * phi2 * fakeError) + dePhredTable[phi1]);
}

void TSimulatorReadBQSR::initializeFakeQualToTrueQualTable(){
	if(!fakeQualToTrueQualTableInitialized){
		fakeQualToTrueQual = new double[maxQualPlusOne];
		for(int i=0; i<maxQualPlusOne; ++i)
			fakeQualToTrueQual[i] = returnTrueQual(i);
	}
	fakeQualToTrueQualTableInitialized = true;
};

void TSimulatorReadBQSR::calculateSlopeIntercept(){
	double sum = 0.0;
	double normalization_sum = 0.0;
	//gamma density starts at 0 but p at 1!

	for(int p=1; p<(readLengthDist->max() + 1) ; ++p){
		sum += (double) p * (1 - readLengthDist->gammaCumulDensity[p-1]);
//		std::cout << "p " << p << " (1 - readLengthDist->gammaCumulDensity[p-1]) / normalization_sum " << (1 - readLengthDist->gammaCumulDensity[p-1]) / normalization_sum << std::endl;
	}
	m = (1.0 - revIntercept) / sum;
	intercept = revIntercept - m * readLengthDist->max();
}

double TSimulatorReadBQSR::returnBetaPp(int pos){
	return(m * (double) pos + intercept);
}

void TSimulatorReadBQSR::simulate(short* posAddress, readLengthContainer & rl, TGenotypeMap & genoMap){
	static short base;
	static int fakeQual, trueQual;
	static long p;
	queryBases = "";
	bamQualities = "";

	for(p=0; p<rl.readLength; ++p){
		//get true nucleotide
		base = *(posAddress + p);

		//apply PMD
		if(pmdInitialized){
			applyPMD(base, p, rl);
		}
		//sample quality and add error
		fakeQual = sampleFakeQuality();
		if(fakeQual > maxQual) fakeQual = (char) maxQual;

		trueQual = fakeQualToTrueQual[fakeQual];
		if(randomGenerator->getRand() < (qualToErroTable[trueQual] * returnBetaPp(p))){ //qualToError knows that quals are in ascii
			base = (base + randomGenerator->pickOne(3) + 1) % 4;
		}
		//add to bam alignment
		bamQualities += (char) (fakeQual + 33);
		queryBases += toBase[base];
	}
};

//--------------------------
//BQSR quality transformation
//-------------------------
TSimulatorReadBQSRQual::TSimulatorReadBQSRQual(TSimulatorReadLength* ReadLengthDist, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase): TSimulatorReadBQSR(ReadLengthDist, params, Logfile, RandomGenerator, ToBase){
	//all beta_p should be 1
	intercept = 1.0;
	revIntercept = 1.0;
	m = 1.0;
}

TSimulatorReadBQSRPos::TSimulatorReadBQSRPos(TSimulatorReadLength* ReadLengthDist, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase): TSimulatorReadBQSR(ReadLengthDist, params, Logfile, RandomGenerator, ToBase){
	revIntercept = params.getParameterDoubleWithDefault("BQSRPosition", 2.0);
	if(revIntercept < 0) throw("BQSRPosition cannot be negative!");
	logfile->list("ReverseIntercept is set to " + toString(revIntercept));
	calculateSlopeIntercept();
	logfile->conclude("Simulating BQSR position effect of quality distortion with slope  = " + toString(m) + " and intercept = " + toString(intercept));
	if(intercept < 0) throw("Intercept for position transformation is negative -> choose higher value for BQSRPosition!");
}









