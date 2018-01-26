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
	size_t pos = s.find("(");
	if(pos == std::string::npos)
		_max = stringToIntCheck(s);
	else if(pos == 0){
		pos = s.find(')');
		if(pos == std::string::npos || pos != s.size() - 1)
			throw "Failed to understand fixed quality '" + s + "'!";
		_max = stringToIntCheck(s.substr(1,pos - 1));
	} else
		throw "Failed to understand fixed quality '" + s + "'!";


	_min = -1;
	tmpInt = 0;
	_maxPlusOne = -1;
	_mean = -1.0;
	_sd = -1.0;
};
TSimulatorQualityDist::TSimulatorQualityDist(){
	_max = 30;
	_min = -1;
	tmpInt = 0;
	_maxPlusOne = -1;
	_mean = -1.0;
	_sd = -1.0;
}

void TSimulatorQualityDist::sample(int* qualities, const int & len){
	for(tmpInt=0; tmpInt<len; ++tmpInt)
		qualities[tmpInt] = _max;
};

void TSimulatorQualityDist::printDetails(TLog* logfile){
	logfile->list("Fixed quality of " + toString(_max));
};

TSimulatorQualityDistBinned::TSimulatorQualityDistBinned(std::string & s, TRandomGenerator* RandomGenerator):TSimulatorQualityDist(){
	size_t pos = s.find("(");
	if(pos == 0){
		s.erase(0,1);
		pos = s.find(')');
		if(pos == std::string::npos || pos != s.size() - 1)
			throw "Failed to understand binned quality '" + s + "'! Use binned(quality_1,quality_2,..,quality_n).";
		s.erase(pos,1);
		fillVectorFromString(s, qualBins, ',');
	} else
		throw "Failed to understand binned quality '" + s + "'! Use binned(quality_1,quality_2,..,quality_n).";

	numQualBins = qualBins.size();
	randomGenerator = RandomGenerator;

}

void TSimulatorQualityDistBinned::sample(int* qualities, const int & len){
	for(int i=0; i<len; ++i){
		qualities[i] = qualBins[randomGenerator->pickOne(numQualBins)];
	}
}

void TSimulatorQualityDistBinned::printDetails(TLog* logfile){
	std::string tmpS = concatenateString(qualBins, ",");
	logfile->list("Quality scores uniformly distributed among the following bins: " + tmpS);
};

TSimulatorQualityDistNormal::TSimulatorQualityDistNormal(std::string & s, TRandomGenerator* RandomGenerator):TSimulatorQualityDist(){
	densitiesInitialized = false;
	randomGenerator = RandomGenerator;

	parseFunctionString(s);
	_maxPlusOne = _max + 1;

	fillDensities();
};

TSimulatorQualityDistNormal::TSimulatorQualityDistNormal(double & mean, double & sd, int min, int max, TRandomGenerator* RandomGenerator):TSimulatorQualityDist(){
	densitiesInitialized = false;
	randomGenerator = RandomGenerator;

	_mean = mean;
	_sd = sd;
	_min = min;
	_max = max;
	_maxPlusOne = _max + 1;

	fillDensities();
};


void TSimulatorQualityDistNormal::parseFunctionString(std::string & s){
	std::string orig = s;

	if(s[0] != '(')
		throw "Fail to understand function '" + orig + "': use format normal(mean,sd)[min,max].";
	s.erase(0,1);

	unsigned int pos = s.find(",");
	if(pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format normal(mean,sd)[min,max].";
	_mean = stringToDouble(s.substr(0,pos));
	if(_mean < 0)
		throw "Fail to understand function '" + orig + "': mean must be > 0.";
	s.erase(0,pos+1);

	pos = s.find(")");
	if(pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format normal(mean,sd)[min,max].";
	_sd = stringToDouble(s.substr(0,pos));
	if(_sd < 0)
			throw "Fail to understand function '" + orig + "': sd must be > 0.";
	s.erase(0,pos+1);

	if(s[0] != '[')
		throw "Fail to understand function '" + orig + "': use format normal(mean,sd)[min,max].";
	s.erase(0,1);
	pos = s.find(",");
	if(pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format normal(mean,sd)[min,max].";
	_min = stringToDouble(s.substr(0,pos));
	if(_min < 0)
		throw "Fail to understand function '" + orig + "': min must be >= 0!";
	s.erase(0,pos+1);
	pos = s.find("]");
	if(pos == std::string::npos)
		throw "Fail to understand function '" + orig + "': use format normal(mean,sd)[min,max].";
	_max = stringToDouble(s.substr(0,pos));
	if(_max < _min)
			throw "Fail to understand function '" + orig + "': max must be >= min!";
};

void TSimulatorQualityDistNormal::fillDensities(){
	//fill densities
	size = _maxPlusOne - _min;
	densities = new double[size];
	cumulDensities = new double[size];

	double nextDens = randomGenerator->normalCumulativeDistributionFunction(_min-0.5, _mean, _sd);
	double prevDens;
	double sum = 0;
	for(int i=0; i<size; ++i){
		prevDens = nextDens;
		nextDens = randomGenerator->normalCumulativeDistributionFunction(_min + i + 0.5, _mean, _sd);
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

void TSimulatorQualityDistNormal::sample(int* qualities, const int & len){
	for(tmpInt=0; tmpInt<len; ++tmpInt){
		qualities[tmpInt] = randomGenerator->pickOne(size, cumulDensities) + _min;
	}
};

void TSimulatorQualityDistNormal::printDetails(TLog* logfile){
	logfile->list("Normally distributed quality scores with mean=" + toString(_mean) + " and sd=" + toString(_sd) + ", truncated to [" + toString(_min) + "," + toString(_max) + "].");
};


//-----------------------------------------------
//TSimulatorQualityTransformation
//-----------------------------------------------
TSimulatorQualityTransformation::TSimulatorQualityTransformation(TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator){
	qualityDist = QualityDist;
	randomGenerator = RandomGenerator;
	p = 0;
};


void TSimulatorQualityTransformation::simulateQualitiesAndErrors(Base* bases, int* qualities, const int & len){
	//simulate qualities
	qualityDist->sample(qualities, len);

	//add errors
	for(p=0; p<len; ++p){
		if(randomGenerator->getRand() < qualityMap.phredIntToErrorMap[qualities[p]])
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
		throw "Wrong number of beta values for recal quality transformation (" + toString(betas.size()) + " instead of 24)! Require one for quality, quality^2, position, position^2 and one each for all 20 contexts.";

	//precalculate stuff
	fillTransformationTable(maxReadLength);
};

void TSimulatorQualityTransformationRecal::fillTransformationTable(int maxReadLength){
	//set size
	maxReadLengthPlusOne = maxReadLength + 1;
	maxQualPlusOne = qualityDist->max() + 1;
	numContext = genoMap.numContextsNotN;

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
			transformedQuality[q][p] = new int[numContext];
			for(c=0; c<numContext; ++c){
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
				transformedQuality[q][p][c] = round(-10.0 * log10(transQual / (1.0 + transQual)));
			}
		}
	}

	//clean up
	delete[] qualTermForTransformation;
	delete[] posTermForTransformation;
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

void TSimulatorQualityTransformationRecal::simulateQualitiesAndErrors(Base* bases, int* qualities, const int & len){
	//simulate qualities
	qualityDist->sample(qualities, len);

	//add errors and transform qualities
	previousBase = N;
	for(p=0; p<len; ++p){
		if(randomGenerator->getRand() < qualityMap.phredIntToErrorMap[qualities[p]])
			bases[p] = static_cast<Base>( (bases[p] + randomGenerator->pickOne(3) + 1) % 4);

		//transform qualities
		qualities[p] = transformedQuality[qualities[p]][p][genoMap.contextMap[previousBase][bases[p]]];
		previousBase = bases[p];
	}
};



//------------------------------------
//TSimulatorQualityTransformationBQSR
//------------------------------------
TSimulatorQualityTransformationBQSR::TSimulatorQualityTransformationBQSR(const std::string & s, TSimulatorReadLength* ReadLengthDist, TLog* logfile, TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator): TSimulatorQualityTransformation(QualityDist, RandomGenerator){
	readLengthDist = ReadLengthDist;
	maxReadLength = readLengthDist->max();
	minPhredInt = 1;
	maxPhredInt = 42;
	maxPhredIntPlusOne = maxPhredInt + 1;
	meanPhred = qualityDist->mean();
	sdPhred = qualityDist->sd();
	std::vector<std::string> vec;
	fillVectorFromStringAnySkipEmpty(s, vec, ",");
	phi1 = stringToInt(vec[0]);
	phi2 = stringToDouble(vec[1]);
	logfile->list("Simulating BQSR quality effect with phi1 = " + toString(phi1) + " and phi2 = " + toString(phi2));

	revIntercept = stringToDouble(vec[2]);
	if(revIntercept < 0.0) throw("revIntercept cannot be negative!");
	else if(revIntercept == 1.0){
		m = 0.0;
		intercept = 1.0;
		logfile->list("revIntercept is set to 1 -> Simulating only BQSR quality effect.");
	} else {
		calculateSlopeIntercept();

		//check if beta_q * beta_p ever results in error > 1 (if yes increase minPhredInt)
		double betaQq = returnTrueError(minPhredInt);
		while (betaQq * revIntercept > 1.0){
			++minPhredInt;
			betaQq = returnTrueError(minPhredInt);
		}

		//report
		logfile->list("Simulating BQSR fake qualities to range from " + toString(minPhredInt) + " to " + toString(maxPhredInt));
		logfile->list("Simulating BQSR position effect of quality distortion with slope  = " + toString(m) + " and intercept = " + toString(intercept));
	}

//	parseBQSRQualInput(params);
	setFakePhredDistribution(logfile); //first find kappa and lambda
	fakePhredDist = new TSimulatorQualityDistNormal(kappa, lambda, minPhredInt, maxPhredInt, RandomGenerator);

}

void TSimulatorQualityTransformationBQSR::calculateSlopeIntercept(){
	double sum = 0.0;
	//gamma density starts at 0 but p at 1!
	for(int p=0; p<maxReadLength ; ++p)
		sum += (double) (p+1) * readLengthDist->positionProbs[p];

	m = (1.0 - revIntercept) / (sum - maxReadLength);
	intercept = revIntercept - m * (double) maxReadLength;

	if(intercept < 0) throw "The value given for the reverse intercept results in a negative intercept!";
}

int TSimulatorQualityTransformationBQSR::sampleFakePhredInt(){
	int qualPhredInt = round(randomGenerator->getNormalRandom(kappa, lambda));
	if(qualPhredInt > 93) qualPhredInt = 93;
	if(qualPhredInt < 0) qualPhredInt = 0;
	return qualPhredInt;
}


void TSimulatorQualityTransformationBQSR::fillQBetaQBetaP(){
	double init_value = -1.0;
	double betaQq;

	QBetaQBetaP.resize( maxPhredIntPlusOne , std::vector<double>( maxReadLength , init_value ) );
	errorBetaQBetaP.resize( maxPhredIntPlusOne , std::vector<double>( maxReadLength , init_value ) );

	for(int q = 0; q < maxPhredIntPlusOne; ++ q){
		betaQq = returnTrueError(q);
		for(int p = 0; p<maxReadLength; ++p){
			errorBetaQBetaP[q][p] = (betaQq) * returnBetaPp(p);
			QBetaQBetaP[q][p] = qualityMap.errorToPhred(errorBetaQBetaP[q][p]);
		}
	}
}

void TSimulatorQualityTransformationBQSR::fillWeights(double & kappa_cur, double & lambda_cur){
	w = new double[maxPhredIntPlusOne];

	//w at minQual
	w[minPhredInt] = randomGenerator->normalCumulativeDistributionFunction(((double) minPhredInt + 0.5), kappa_cur, lambda_cur);

	//w at intermediate Q
	for(int q = (minPhredInt + 1); q < maxPhredInt; ++q){
		double start = randomGenerator->normalCumulativeDistributionFunction((double) q - 0.5, kappa_cur, lambda_cur);
		double end = randomGenerator->normalCumulativeDistributionFunction((double) q + 0.5, kappa_cur, lambda_cur);
		w[q] = end - start;
	}

	//w at maxQual
	w[maxPhredInt] = 1.0 - randomGenerator->normalCumulativeDistributionFunction(((double) maxPhredInt - 0.5), kappa_cur, lambda_cur);
	weightsInitialized = true;
}

double TSimulatorQualityTransformationBQSR::returnCurMean(){
	double curMean = 0.0;

	for(int q = minPhredInt; q < maxPhredIntPlusOne; ++ q){
		double sumP = 0.0;
		for(int p = 0; p<maxReadLength; ++p){
			sumP += QBetaQBetaP[q][p] * readLengthDist->positionProbs[p];
		}
		curMean += sumP * w[q];
	}

	return(curMean);
}

double TSimulatorQualityTransformationBQSR::returnCurSD(double & kappa){
	float lambda = 0.0;
	for(int q = minPhredInt; q < maxPhredIntPlusOne; ++ q){
		for(int p = 0; p<readLengthDist->max(); ++p){
			lambda += ((QBetaQBetaP[q][p] - kappa) * (QBetaQBetaP[q][p] - kappa) * readLengthDist->positionProbs[p] * w[q]);
		}
	}
	return(sqrt(lambda));
}

double TSimulatorQualityTransformationBQSR::returnDelta(double & curMean, double & curSD){
	double delta;
	delta = (curMean - meanPhred)*(curMean - meanPhred) + (curSD - sdPhred)*(curSD - sdPhred);
	return(delta);
}

//---------------------------------

void TSimulatorQualityTransformationBQSR::setFakePhredDistribution(TLog* logfile){
	double kappa_cur = meanPhred;
	double lambda_cur = sdPhred;
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
	out << kappa_cur << "\t" << lambda_cur << "\t" << delta_cur << std::endl;

	int nTurns = 10;
	int nIter = 50;

	logfile->startIndent("Estimating mean (kappa) and sd (lambda) for distorted quality score distribution");

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


	//	logfile->list("Current estimates: kappa = " + toString(kappa_cur) + ", lambda = " + toString(lambda_cur) + ", delta = " + toString(delta_cur) + " corresponding to a true qual dist = N(" + toString(mean_cur) + "," + toString(sd_cur) + ")");
		out << kappa_cur << "\t" << lambda_cur << "\t" << delta_cur << std::endl;

	}

	logfile->endIndent();

	logfile->conclude("The final estimates for kappa = " + toString(kappa_cur) + " and lambda = " + toString(lambda_cur) + ". This result in a true quality score distribution = N(" + toString(mean_cur) + "," + toString(sd_cur) + ") with delta = " + toString(delta_cur) + ".");
	if(delta_cur >= 0.25) logfile->warning("Current parameter values for phi1, meanQual and sdQual do not allow for accurate estimation of kappa and lambda!");
	kappa = kappa_cur;
	lambda = lambda_cur;
}

double TSimulatorQualityTransformationBQSR::returnTrueError(const int & fakePhredInt){
	return(pow(10.0, -1.0/10.0 * phi2 * (double) fakePhredInt) + qualityMap.phredIntToErrorMap[phi1]);
}

double TSimulatorQualityTransformationBQSR::returnBetaPp(const int & pos){
	return(m * (double) (pos+1) + intercept);
}

void TSimulatorQualityTransformationBQSR::simulateQualitiesAndErrors(Base* bases, int* phredIntQualities, const int & len){
	//for loop that simulates errors according to true qual and returns the fake qualities for bam file
	fakePhredDist->sample(phredIntQualities,len);

	for(p=0; p<len; ++p){
		if(randomGenerator->getRand() < errorBetaQBetaP[phredIntQualities[p]][p]){
			bases[p] = static_cast<Base>( (bases[p] + randomGenerator->pickOne(3) + 1) % 4);
		}
	}
}





