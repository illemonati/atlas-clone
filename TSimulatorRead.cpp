/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#include "TSimulatorRead.h"
#include "TParameters.h"

double TSimulatorRead::dePhredAscii(double x){
	return pow(10, -(x-33.0) / 10.0);
}

double TSimulatorRead::dePhred(double x){
	return pow(10, -(x) / 10.0);
}

int TSimulatorRead::phred(double x){
	return(round(-10.0 * log10(x)));
}

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
	double mQ = params.getParameterDoubleWithDefault("meanQual", 30);
	double sdQ = params.getParameterDoubleWithDefault("sdQual", 10);
	int maxQual = params.getParameterDoubleWithDefault("maxQual", 127);

	setTrueQualityDistribution(mQ, sdQ);
	logfile->list("Will simulate normal distributed quality scores with mean = " + toString(meanQual) + " and sd = " + toString(sdQual) + ", capped at " + toString(maxQual) + ".");

	initializeQualToErrorTable();

};

std::string TSimulatorRead::getQueryBases(){
	return(queryBases);
}
std::string TSimulatorRead::getBamQualities(){
	return(bamQualities);
}

void TSimulatorRead::setTrueQualityDistribution(double mean, double sd){
	meanQual = mean + 33.0; //add 33 to mean quality to get in in char
	sdQual = sd;
}

int TSimulatorRead::sampleTrueQuality(){
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
		for(int i=33; i<maxQual+33; ++i)
			qualToErroTable[i] = dePhredAscii(i);
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
		else bamQualities += (char) transQual;
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
	for(int i=0; i<34; ++i)
		qualTermForTransformation[i] = 100.0;

	double tmp;
	for(int i=34; i<127; ++i){
		tmp = pow(10.0, -(double) (i - 33.0) / 10.0);
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
	return -10.0 * log10(tmp / (1.0 + tmp)) + 33.0;
}

int TSimulatorReadRecal::returnQual(int & qual, int & pos, int & baseContext){
	return transformQuality(qual, pos, baseContext);
}


//--------------------------
//BQSR transformation Position
//-------------------------

TSimulatorReadBQSRPos::TSimulatorReadBQSRPos(TSimulatorReadLength* ReadLengthDist, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase): TSimulatorRead(params, Logfile, RandomGenerator, ToBase){
	//position parameters
	readLengthDist = ReadLengthDist;
	revIntercept = params.getParameterDoubleWithDefault("BQSRPosition", 2.0);
	logfile->list("ReverseIntercept is set to " + toString(revIntercept));
	calculateSlopeIntercept();
	logfile->conclude("Simulating BQSR position effect of quality distortion with slope  = " + toString(m) + " and intercept = " + toString(intercept));

	//quality parameters
	parseBQSRQualInput(params);
	setFakeQualityDistribution(); //first find kappa and lambda
	initializeFakeQualToTrueQualTable();

}

void TSimulatorReadBQSRPos::parseBQSRQualInput(TParameters & params){
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

int TSimulatorReadBQSRPos::sampleFakeQuality(){
	int qual = round(randomGenerator->getNormalRandom(kappa, lambda));
	if(qual > 126) qual = 126;
	if(qual < 33) qual = 33;
	return qual;
}

void TSimulatorReadBQSRPos::setFakeQualityDistribution(){
	double kappa = -1;
	double lambda = -1;
}

int TSimulatorReadBQSRPos::returnTrueQual(int & fakeQual){
	double fakeError = dePhred(fakeQual);
	return(pow(10, -1/10 * phi2 * fakeError) + phred(phi1));
}

void TSimulatorReadBQSRPos::initializeFakeQualToTrueQualTable(){
	if(!fakeQualToTrueQualTableInitialized){
		fakeQualToTrueQual = new double[127];
		for(int i=0; i<33; ++i)
			fakeQualToTrueQual[i] = 1.0;
		for(int i=33; i<maxQual+33; ++i)
			fakeQualToTrueQual[i] = returnTrueQual(i);
	}
	fakeQualToTrueQualTableInitialized = true;
};

void TSimulatorReadBQSRPos::calculateSlopeIntercept(){
	double sum = 0.0;
	for(int p=0; p<readLengthDist->max(); ++p){
		sum += (double) p * readLengthDist->gammaCumulDensity[p];
	}
	m = (1.0 - revIntercept) / sum;
	intercept = revIntercept - m * readLengthDist->max();
}

double TSimulatorReadBQSRPos::returnBetaPp(int pos){
	return(m * (double) pos + intercept);
}

void TSimulatorReadBQSRPos::simulate(short* posAddress, readLengthContainer & rl, TGenotypeMap & genoMap){
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
		if(fakeQual > 127)	fakeQual = (char) maxQual;

		trueQual = fakeQualToTrueQual[fakeQual];
		if(randomGenerator->getRand() < (qualToErroTable[trueQual] * returnBetaPp(p))){ //qualToError knows that quals are in ascii
			base = (base + randomGenerator->pickOne(3) + 1) % 4;
		}
		//add to bam alignment
		bamQualities += (char) fakeQual;
		queryBases += toBase[base];
	}
};

//--------------------------
//BQSR transformation
//-------------------------
TSimulatorReadBQSRQual::TSimulatorReadBQSRQual(TSimulatorReadLength* ReadLengthDist, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase): TSimulatorReadBQSRPos(ReadLengthDist, params, Logfile, RandomGenerator, ToBase){
	//set position parameters to 0
	intercept = 0.0;
	revIntercept = 0.0;
	m = 1.0;
}

/*
void TSimulatorReadBQSRQual::simulate(short* posAddress, readLengthContainer & rl, TGenotypeMap & genoMap){
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
		if(fakeQual > 127)	fakeQual = (char) maxQual;

		trueQual = fakeQualToTrueQual[fakeQual];
		if(randomGenerator->getRand() < qualToErroTable[trueQual]){ //qualToError knows that quals are in ascii
			base = (base + randomGenerator->pickOne(3) + 1) % 4;
		}
		//add to bam alignment
		bamQualities += (char) fakeQual;
		queryBases += toBase[base];
	}
}



*/















