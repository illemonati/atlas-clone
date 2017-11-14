/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#include "TSimulatorRead.h"
#include "TParameters.h"


double TSimulatorRead::dePhred(int x){
	return(pow(10, -((double) x) / 10.0));
}

int TSimulatorRead::phred(double x){
	return(round(-10.0 * log10(x)));
}

void TSimulatorRead::initializeDePhredTable(){
	if(!dePhredTableInitialized){
		dePhredTable = new double[maxQualPlusOne-minQual]; //-minQual
		for(int i=minQual; i<maxQualPlusOne; ++i){
			dePhredTable[i] = dePhred(i);
		}
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
			qualToErroTable[i] = dePhredTable[i];
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
	if(params.getParameterString("BQSRPosition", false) != ""){
		revIntercept = params.getParameterDoubleWithDefault("BQSRPosition", 2.0);
		std::cout << "revIntercept " << revIntercept << std::endl;
		if(revIntercept < 0) throw("BQSRPosition cannot be negative!");
		logfile->list("ReverseIntercept is set to " + toString(revIntercept));
		calculateSlopeIntercept();
		logfile->conclude("Simulating BQSR position effect of quality distortion with slope  = " + toString(m) + " and intercept = " + toString(intercept));
		if(intercept < 0) throw("Intercept for position transformation is negative -> choose higher value for BQSRPosition!");
	} else logfile->list("BQSRPosition parameter not found. Simulating only BQSR quality effect.");

	//quality parameters
	parseBQSRQualInput(params);
	setFakeQualityDistribution(); //first find kappa and lambda
	initializeDePhredTable();
	initializeTrueQualToFakeQualTable();
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

	logfile->conclude("Will simulate quality distortion with alpha1 = " + toString(phi1) + " and alpha2 = " + toString(phi2));
}


void TSimulatorReadBQSR::calculateSlopeIntercept(){
	double sum = 0.0;
	//gamma density starts at 0 but p at 1!
	for(int p=1; p<(readLengthDist->max() + 1) ; ++p){
		sum += (double) p * (1 - readLengthDist->gammaCumulDensity[p-1]);
	}
	m = (1.0 - revIntercept) / sum;
	intercept = revIntercept - m * readLengthDist->max();
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
	double betaQq;

	QBetaQBetaP.resize( num_of_col , std::vector<double>( num_of_row , init_value ) );
	for(int q = minQual; q < maxQualPlusOne; ++ q){
		betaQq = returnFakeError(q);
		for(int p = 0; p<readLengthDist->max(); ++p){
			QBetaQBetaP[q][p] = phred(dePhredTable[phred(betaQq)] * returnBetaPp(p));
		}
	}
	betaQBetaPInitialized = true;
}

void TSimulatorReadBQSR::fillWeights(double & kappa_cur, double & lambda_cur){
	w = new double[maxQualPlusOne - minQual];
	//w at minQual
	w[0] = randomGenerator->normalCumulativeDistributionFunction(((double) minQual + 0.5), kappa_cur, lambda_cur);

	//w at intermediate Q
	for(int q = (minQual + 1); q < maxQual; ++q){
		double start = randomGenerator->normalCumulativeDistributionFunction((double) q - 0.5, kappa_cur, lambda_cur);
		double end = randomGenerator->normalCumulativeDistributionFunction((double) q + 0.5, kappa_cur, lambda_cur);
		w[q-minQual] = end - start;
	}

	//w at maxQual
	w[(maxQual - minQual)] = 1 - randomGenerator->normalCumulativeDistributionFunction(((double) maxQual - 0.5), kappa_cur, lambda_cur);
	weightsInitialized = true;
}

double TSimulatorReadBQSR::returnCurMean(){
	double curMean = 0.0;

	for(int q = minQual; q < maxQualPlusOne; ++ q){
		double sumP = 0.0;
		for(int p = 0; p<readLengthDist->max(); ++p){
			sumP += QBetaQBetaP[q][p] * readLengthDist->positionProbs[p];
		}
		curMean += sumP * w[q];
//		if(q == maxQual){
//			std::cout << "curMean " << curMean << std::endl;
//			std::cout << "sumP " << sumP << std::endl;
//			std::cout << "w[q] " << w[q] << std::endl;
//		}
	}

	return(curMean);
}

double TSimulatorReadBQSR::returnCurSD(double & kappa){
	float lambda = 0.0;
	for(int q = minQual; q < maxQualPlusOne; ++ q){
		for(int p = 0; p<readLengthDist->max(); ++p){
			lambda += ((QBetaQBetaP[q][p] - kappa) * (QBetaQBetaP[q][p] - kappa) * readLengthDist->positionProbs[p] * w[q]);
		}
	}
	return(sqrt(lambda));
}

double TSimulatorReadBQSR::returnDelta(double & curMean, double & curSD){
	double delta;
	delta = (curMean - meanQual)*(curMean - meanQual) + (curSD - sdQual)*(curSD - sdQual);
	std::cout << "delta in returnDelta: " << delta << std::endl;
	return(delta);
}

//---------------------------------

void TSimulatorReadBQSR::setFakeQualityDistribution(){
	double kappa_cur = meanQual;
	double lambda_cur = sdQual;
	double delta_cur, delta_old;
	double mean_cur, sd_cur;
	float stepSize;

	std::ofstream out;
	std::string filename = "path_simulation.txt";
	out.open(filename.c_str());
	out << "kappa\tlambda\tdelta\n";

	fillQBetaQBetaP();
	fillWeights(kappa_cur, lambda_cur);
	mean_cur = returnCurMean();
	sd_cur = returnCurSD(mean_cur);

	delta_cur = returnDelta(mean_cur, sd_cur);
	std::cout << "initial values: kappa lambda error:  "<< kappa_cur << "\t" << lambda_cur << "\t" << delta_cur << std::endl;

	out << kappa_cur << "\t" << lambda_cur << "\t" << delta_cur << std::endl;

	int nTurns = 10;
	int nIter = 50;

	logfile->startIndent("Estimating mean (kappa) and sd (lambda) for distorted quality score distribution:");

	//optimize one param at a time, update, optimize again

	for(int t=0; t<nTurns; ++t){
		std::cout << "####### new turn" << std::endl;
		//update kappa
		stepSize = 5.0;
		for(int i=0; i<nIter; ++i){
			std::cout << "####### new iter kappa" << std::endl;

			//move and calc error
			delta_old = delta_cur;
			kappa_cur += stepSize;
			fillWeights(kappa_cur, lambda_cur);
			mean_cur = returnCurMean();
			sd_cur = returnCurSD(mean_cur);
			delta_cur = returnDelta(mean_cur, sd_cur);
			if(delta_cur > delta_old)
				stepSize = -stepSize/exp(1);
		}

		//update lambda
		stepSize = 1.0;
		for(int i=0; i<nIter; ++i){
			std::cout << "####### new iter lambda" << std::endl;

			//move and calc error
			delta_old = delta_cur;
			lambda_cur += stepSize;
			fillWeights(kappa_cur, lambda_cur);
			mean_cur = returnCurMean();
			sd_cur = returnCurSD(mean_cur);
			delta_cur = returnDelta(kappa_cur, sd_cur);
			if(delta_cur > delta_old)
				stepSize = -stepSize/exp(1);
		}


		logfile->list("Current estimates: kappa = " + toString(kappa_cur) + ", lambda = " + toString(lambda_cur) + ", delta = " + toString(delta_cur));
		out << kappa_cur << "\t" << lambda_cur << "\t" << delta_cur << std::endl;
	}

	logfile->conclude("The final estimates for kappa and lambda result in a true quality score being simulated according to N(" + toString(returnCurMean()) + "," + toString(returnCurSD(kappa_cur)) + "). This corresponds to a delta of " + toString(delta_cur) + ".");
	if(delta_cur >= 0.25) logfile->warning("Current parameter values for phi1, meanQual and sdQual do not allow for accurate estimation of kappa and lambda!");
	kappa = kappa_cur;
	lambda = lambda_cur;

	logfile->endIndent();
}

double TSimulatorReadBQSR::returnFakeError(int & trueQual){
	return(pow(10, -1/10.0 * phi2 * trueQual) + dePhredTable[phi1]);
}

void TSimulatorReadBQSR::initializeTrueQualToFakeQualTable(){
	if(!trueQualToFakeQualTableInitialized){
		trueQualToFakeQual = new double[maxQualPlusOne];
		for(int i=0; i<maxQualPlusOne; ++i){
			trueQualToFakeQual[i] = round(phred(returnFakeError(i)));
		}
	}
	trueQualToFakeQualTableInitialized = true;
};



double TSimulatorReadBQSR::returnBetaPp(int & pos){
	return(m * (double) pos + intercept);
}

void TSimulatorReadBQSR::simulate(short* posAddress, readLengthContainer & rl, TGenotypeMap & genoMap){
	static short base;
	static int fakeQual, QBetaQq;
	static long p;
	int pInt;
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

		QBetaQq = trueQualToFakeQual[fakeQual];
		double BetaQp;
		pInt = (int) p;
		BetaQp = returnBetaPp(pInt);
		if(randomGenerator->getRand() < (qualToErroTable[QBetaQq] * BetaQp)){ //qualToError knows that quals are in ascii
			base = (base + randomGenerator->pickOne(3) + 1) % 4;
		}
		//add to bam alignment
		bamQualities += (char) (fakeQual + 33);
		queryBases += toBase[base];
	}
};









