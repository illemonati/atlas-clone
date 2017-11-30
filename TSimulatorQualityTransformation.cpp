/*
 * TSimulatorQualityTransformation.cpp
 *
 *  Created on: Nov 28, 2017
 *      Author: phaentu
 */


#include "TSimulatorQualityTransformation.h"


//----------------------------------
//TSimulatorQualityDist
//----------------------------------
TSimulatorQualityDist::TSimulatorQualityDist(std::string & s){
	_max = stringToIntCheck(s);
	tmpInt = 0;
};
TSimulatorQualityDist::TSimulatorQualityDist(){
	_max = 30;
	tmpInt = 0;
}

void TSimulatorQualityDist::sample(int* qualities, int & len){
	for(tmpInt=0; tmpInt<len; ++tmpInt)
		qualities[tmpInt] = _max;
};

void TSimulatorQualityDist::printDetails(TLog* logfile){
	logfile->list("Fixed quality of " + toString(_max));
};

TSimulatorQualityDistNormal::TSimulatorQualityDistNormal(std::string & s, TRandomGenerator* RandomGenerator):TSimulatorQualityDist(){
	densitiesInitialized = false;
	randomGenerator = RandomGenerator;

	parseFunctionString(s);
	_maxPlusOne = _max + 1;

	fillDensities();
};

void TSimulatorQualityDistNormal::parseFunctionString(std::string & s){
	std::string orig = s;

	if(s[0] != '(')
		throw "Fail to understand function '" + orig + "': use format normal(var1,var2)[min,max].";
	s.erase(0,1);

	unsigned int pos = s.find(",");
	if(pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format normal(var1,var2)[min,max].";
	mean = stringToDouble(s.substr(0,pos));
	if(mean < 0)
		throw "Fail to understand function '" + orig + "': mean must be > 0.";
	s.erase(0,pos+1);

	pos = s.find(")");
	if(pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format normal(var1,var2)[min,max].";
	sd = stringToDouble(s.substr(0,pos));
	if(sd < 0)
			throw "Fail to understand function '" + orig + "': sd must be > 0.";
	s.erase(0,pos+1);

	if(s[0] != '[')
		throw "Fail to understand function '" + orig + "': use format normal(var1,var2)[min,max].";
	s.erase(0,1);
	pos = s.find(",");
	if(pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format normal(var1,var2)[min,max].";
	_min = stringToDouble(s.substr(0,pos));
	if(_min < 0)
		throw "Fail to understand function '" + orig + "': min must be >= 0!";
	s.erase(0,pos+1);
	pos = s.find("]");
	if(pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format normal(var1,var2)[min,max].";
	_max = stringToDouble(s.substr(0,pos));
	if(_max < _min)
			throw "Fail to understand function '" + orig + "': max must be >= min!";
};

void TSimulatorQualityDistNormal::fillDensities(){
	//fill densities
	size = _maxPlusOne - _min;
	densities = new double[size];
	cumulDensities = new double[size];

	double nextDens = randomGenerator->normalCumulativeDistributionFunction(_min-0.5, mean, sd);
	double prevDens;
	double sum = 0;
	for(int i=0; i<size; ++i){
		prevDens = nextDens;
		nextDens = randomGenerator->normalCumulativeDistributionFunction(_min + i + 0.5, mean, sd);
		densities[i] =  nextDens - prevDens;
		sum += densities[i];
	}

	//normalize
	for(int i=0; i<size; ++i)
		densities[i] /= sum;

	//fill cumulative
	cumulDensities[0] = densities[0];
	for(int i=1; i<size; ++i)
		cumulDensities[i] = cumulDensities[i-1] + densities[i];
	cumulDensities[size-1] = 1.0;

	densitiesInitialized = true;
}

int TSimulatorQualityDistNormal::sample(){
	tmpInt = randomGenerator->pickOne(size, cumulDensities);
	return tmpInt + _min;
};

void TSimulatorQualityDistNormal::sample(int* qualities, int & len){
	for(tmpInt=0; tmpInt<len; ++tmpInt){
		qualities[tmpInt] = randomGenerator->pickOne(size, cumulDensities) + _min;
	}
};

void TSimulatorQualityDistNormal::printDetails(TLog* logfile){
	logfile->list("Normally distributed quality scores with mean=" + toString(mean) + " and sd=" + toString(sd) + ", truncated to [" + toString(_min) + "," + toString(_max) + "].");
};


//-----------------------------------------------
//TSimulatorQualityTransformation
//-----------------------------------------------
TSimulatorQualityTransformation::TSimulatorQualityTransformation(TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator){
	qualityDist = QualityDist;
	randomGenerator = RandomGenerator;
	p = 0;
};


void TSimulatorQualityTransformation::simulateQualitiesAndErrors(Base* bases, int* qualities, int & len){
	//simulate qualities
	qualityDist->sample(qualities, len);

	//add errors
	for(p=0; p<len; ++p){
		if(randomGenerator->getRand() < qualityMap.phredToErrorMap[qualities[p]])
			bases[p] = static_cast<Base>((bases[p] + randomGenerator->pickOne(3) + 1) % 4);
	}
};

void TSimulatorQualityTransformation::printDetails(TLog* logfile){
	logfile->list("Will write original qualities to BAM file (no transformation).");
};


//-------------------------------------
//TSimulatorQualityTransformationRecal
//-------------------------------------
TSimulatorQualityTransformationRecal::TSimulatorQualityTransformationRecal(const std::string & s, int maxReadLength, TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator)
	:TSimulatorQualityTransformation(QualityDist, RandomGenerator){

	//string of betas, comma seperated and potentially with repeated indexes
	std::vector<std::string> vec;
	fillVectorFromStringAnySkipEmpty(s, vec, ",");
	repeatIndexes(vec, betas);

	//check if 24 betas were provided
	if(betas.size() != 24)
		throw "Wrong number of beta values for quality transformation (" + toString(betas.size()) + " instead of 24)! Require one for quality, quality^2, position, position^2 and one each for all 20 contexts.";

	//precalculate stuff
	fillTransformationTable(maxReadLength);
};

void TSimulatorQualityTransformationRecal::fillTransformationTable(int maxReadLength){
	//set size
	maxReadLengthPlusOne = maxReadLength + 1;
	maxQualPlusOne = qualityDist->max() + 1;
	numContextPlusOne = genoMap.numContexts + 1;

	//quality term
	double* qualTermForTransformation = new double[maxQualPlusOne];
	double tmp;
	for(int i=0; i<maxQualPlusOne; ++i){
		tmp = pow(10.0, -(double) i / 10.0);
		qualTermForTransformation[i] = log(tmp / (1.0 - tmp));
	}

	//position term
	double* posTermForTransformation = new double[maxReadLengthPlusOne];
	for(int i=0; i<maxReadLengthPlusOne; ++i){
		posTermForTransformation[i] = betas[2] * i + betas[3] * i*i;
	}

	//tmp variables
	int q, p, c;
	double constant;
	double transQual;

	//now fill table
	transformedQuality = new int**[maxQualPlusOne];
	for(q=0; q<maxQualPlusOne; ++q){
		transformedQuality[q] = new int*[maxReadLengthPlusOne];
		for(p=0; p<maxReadLengthPlusOne; ++p){
			transformedQuality[q][p] = new int[numContextPlusOne];
			for(c=0; c<numContextPlusOne; ++c){

				//now calc transformed quality
				constant = posTermForTransformation[p] + betas[c+4] - qualTermForTransformation[q];
				if(4.0 * betas[1] * constant > betas[0] * betas[0]) throw "beta[0]^2 cannot be smaller than 4beta[1](position + context constants)";
				if(betas[1] == 0.0){
					transQual = -constant / betas[0];
				} else {
					tmp = sqrt(betas[0] * betas[0] - 4.0 * betas[1] * constant);
					transQual = (tmp - betas[0]) / 2.0 / betas[1];
					//if(q < 0) q = (-tmp - beta[0]) / 2.0 / beta[1];
				}

				transQual = exp(transQual);
				if(transQual == 0) throw "choose different quality transformation parameters! transQual == 0";
				transformedQuality[q][p][c] = -10.0 * log10(transQual / (1.0 + transQual));
			}
		}
	}

};

void TSimulatorQualityTransformationRecal::clearTransformationTable(){
	for(int q=0; q<maxQualPlusOne; ++q){
		for(p=0; p<maxReadLengthPlusOne; ++p)
			delete[] transformedQuality[q][p];
		delete[] transformedQuality[q];
	}
	delete[] transformedQuality;
};

void TSimulatorQualityTransformationRecal::printDetails(TLog* logfile){
	std::string s = concatenateString(betas, ",");
	logfile->list("Will transform qualities using the recal model with beta = [" + s + "]");
};

void TSimulatorQualityTransformationRecal::simulateQualitiesAndErrors(Base* bases, int* qualities, int & len){
	//simulate qualities
	qualityDist->sample(qualities, len);

	//add errors and transform qualities
	previousBase = N;
	for(p=0; p<len; ++p){
		if(randomGenerator->getRand() < qualityMap.phredToErrorMap[qualities[p]])
			bases[p] = static_cast<Base>( (bases[p] + randomGenerator->pickOne(3) + 1) % 4);

		//transform qualities
		qualities[p] = transformedQuality[qualities[p]][p][genoMap.contextMap[previousBase][bases[p]]];
		previousBase = bases[p];
	}
};


/*
//------------------------------------
//TSimulatorQualityTransformationBQSR
//------------------------------------
TSimulatorQualityTransformationBQSR::TSimulatorQualityTransformationBQSR(TSimulatorReadLength* ReadLengthDist, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator, char* ToBase): TSimulatorRead(params, Logfile, RandomGenerator, ToBase){
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
	} else {
		revIntercept = 1.0;
		intercept = 1.0;
		m = 0.0;
		logfile->list("BQSRPosition parameter not found. Simulating only BQSR quality effect.");
	}

	//quality parameters
	parseBQSRQualInput(params);
	setFakeQualityDistribution(); //first find kappa and lambda
	initializeDePhredTable();
	initializeTrueQualToFakeQualTable();
}

void TSimulatorQualityTransformationBQSR::parseBQSRQualInput(TParameters & params){
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


void TSimulatorQualityTransformationBQSR::calculateSlopeIntercept(){
	double sum = 0.0;
	//gamma density starts at 0 but p at 1!
	for(int p=1; p<(readLengthDist->max() + 1) ; ++p){
		sum += (double) p * readLengthDist->positionProbs[p];
	}
	m = (1.0 - revIntercept) / (sum - readLengthDist->max());
	intercept = revIntercept - m * readLengthDist->max();
}

int TSimulatorQualityTransformationBQSR::sampleFakeQuality(){
	int qual = round(randomGenerator->getNormalRandom(kappa, lambda));
	if(qual > 93) qual = 93;
	if(qual < 0) qual = 0;
	return qual;
}

//---------------------------------
// optimization functions
//---------------------------------

void TSimulatorQualityTransformationBQSR::fillQBetaQBetaP(){
	int num_of_row = maxPos;
	int num_of_col = maxQualPlusOne;
	double init_value = -1.0;
	double betaQq;

	QBetaQBetaP.resize( num_of_col , std::vector<double>( num_of_row , init_value ) );

	for(int q = minQual; q < maxQualPlusOne; ++ q){
		betaQq = returnFakeError(q);
		float a = 0;
		for(int p = 0; p<readLengthDist->max(); ++p){
			QBetaQBetaP[q][p] = phred(dePhredTable[phred(betaQq)] * returnBetaPp(p));
			a += returnBetaPp(p);
		}
	}
	betaQBetaPInitialized = true;
}

void TSimulatorQualityTransformationBQSR::fillWeights(double & kappa_cur, double & lambda_cur){
	w = new double[maxQualPlusOne];

	//w at minQual
	w[minQual] = randomGenerator->normalCumulativeDistributionFunction(((double) minQual + 0.5), kappa_cur, lambda_cur);

	//w at intermediate Q
	for(int q = (minQual + 1); q < maxQual; ++q){
		double start = randomGenerator->normalCumulativeDistributionFunction((double) q - 0.5, kappa_cur, lambda_cur);
		double end = randomGenerator->normalCumulativeDistributionFunction((double) q + 0.5, kappa_cur, lambda_cur);
		w[q] = end - start;
	}

	//w at maxQual
	w[maxQual] = 1 - randomGenerator->normalCumulativeDistributionFunction(((double) maxQual - 0.5), kappa_cur, lambda_cur);
	weightsInitialized = true;
}

double TSimulatorQualityTransformationBQSR::returnCurMean(){
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

double TSimulatorQualityTransformationBQSR::returnCurSD(double & kappa){
	float lambda = 0.0;
	for(int q = minQual; q < maxQualPlusOne; ++ q){
		for(int p = 0; p<readLengthDist->max(); ++p){
			lambda += ((QBetaQBetaP[q][p] - kappa) * (QBetaQBetaP[q][p] - kappa) * readLengthDist->positionProbs[p] * w[q]);
		}
	}
	return(sqrt(lambda));
}

double TSimulatorQualityTransformationBQSR::returnDelta(double & curMean, double & curSD){
	double delta;
	delta = (curMean - meanQual)*(curMean - meanQual) + (curSD - sdQual)*(curSD - sdQual);
	return(delta);
}

//---------------------------------

void TSimulatorQualityTransformationBQSR::setFakeQualityDistribution(){
	double kappa_cur = meanQual;
	double lambda_cur = sdQual;
	double delta_cur, delta_old;
	double mean_cur, sd_cur;
	double stepSize;

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
		//update kappa
		stepSize = 5.0;
		for(int i=0; i<nIter; ++i){
			//move and calc error
			delta_old = delta_cur;
			kappa_cur += stepSize;
			fillWeights(kappa_cur, lambda_cur);
			mean_cur = returnCurMean();
			sd_cur = returnCurSD(mean_cur);
			delta_cur = returnDelta(mean_cur, sd_cur);
			if(delta_cur > delta_old)
				stepSize = -stepSize/exp(1.0);
		}

		//update lambda
		stepSize = 1.0;
		for(int i=0; i<nIter; ++i){
			//move and calc error
			delta_old = delta_cur;
			lambda_cur += stepSize;
			fillWeights(kappa_cur, lambda_cur);
			mean_cur = returnCurMean();
			sd_cur = returnCurSD(mean_cur);
			delta_cur = returnDelta(mean_cur, sd_cur);
			if(delta_cur > delta_old)
				stepSize = -stepSize/exp(1.0);
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

double TSimulatorQualityTransformationBQSR::returnFakeError(int & trueQual){
	return(pow(10, -1/10.0 * phi2 * trueQual) + dePhredTable[phi1]);
}

void TSimulatorQualityTransformationBQSR::initializeTrueQualToFakeQualTable(){
	if(!trueQualToFakeQualTableInitialized){
		trueQualToFakeQual = new double[maxQualPlusOne];
		for(int i=0; i<maxQualPlusOne; ++i){
			trueQualToFakeQual[i] = phred(returnFakeError(i));
		}
	}
	trueQualToFakeQualTableInitialized = true;
};



double TSimulatorQualityTransformationBQSR::returnBetaPp(int & pos){
//	std::cout << "m: " << m << " intercept " << intercept << std::endl;
	return(m * (double) pos + intercept);
}

void TSimulatorQualityTransformationBQSR::simulate(short* posAddress, readLengthContainer & rl, TGenotypeMap & genoMap){
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

*/



