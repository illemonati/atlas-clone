/*
 * SFS.cpp
 *
 *  Created on: Apr 10, 2017
 *      Author: wegmannd
 */

#include "SFS.h"

//--------------------------------
//Class to store and SFS
//--------------------------------
SFS::SFS(int numChr){
	rescaled = false;
	initDimension(numChr + 1);
	monoFrac = 0.0;
}

SFS::SFS(std::string filename){
	storageInitialized = false;
	std::ifstream file(filename.c_str());
	if(!file)
		throw "Failed to open SFS file '" + filename + "'!";
	std::vector<double> vec;
	fillVectorFromLineWhiteSpaceSkipEmptyCheck(file, vec);

	//init dimension
	initDimension(vec.size());

	//now store as fraction
	int i=0;
	float sum = 0;
	for(std::vector<double>::iterator it=vec.begin(); it!=vec.end(); ++it, ++i){
		sfs[i] = *it;
		sum += sfs[i];
	}
	for(i=0; i<dimension; ++i)
		sfs[i] = sfs[i] / sum;
	rescaled = false;
	monoFrac = sfs[0];
	fillCumulative();
}

SFS::SFS(SFS* other, float MonoFrac){
	//construct from other SFS, but rescale SFS to have a specific fraction of monomorphic sites
	storageInitialized = false;
	initDimension(other->dimension);

	//now copy and rescale
	rescaled = true;
	monoFrac = MonoFrac;
	float sum = 0.0;
	sfs[0] = monoFrac;
	for(int i=1; i<dimension; ++i){
		sfs[i] = other->sfs[i];
		sum += sfs[i];
	}
	sum = sum / (1.0 - monoFrac);
	for(int i=1; i<dimension; ++i){
		sfs[i] = sfs[i] / sum;
	}
	fillCumulative();
}

SFS::SFS(int numChr, float theta){
	storageInitialized = false;
	rescaled = false;
	initDimension(numChr  + 1);

	//generate sfs from theta
	float sum = 0;
	sfsFrequencies[0] = 0.0;
	for(int i=1; i<dimension; ++i){
		sfs[i] = theta / (float) i;
		sum += sfs[i];
		sfsFrequencies[i] = (float) i / (float) numChr;
	}
	if(sum > 1.0)
		throw "The choice of theta and sample size results in too many mutations in the SFS!";

	sfs[0] = 1.0 - sum;
	monoFrac = sfs[0];
	fillCumulative();
};

SFS::SFS(int numChr, int onlyThisBin){
	storageInitialized = false;
	rescaled = false;
	initDimension(numChr  + 1);

	//set all to zero except this one bin
	for(int i=0; i<dimension; ++i){
		sfs[i] = 0.0;
	}
	sfs[onlyThisBin] = 1.0;
	monoFrac = sfs[0];
	fillCumulative();
}


void SFS::initDimension(int size){
	dimension = size;
	dimensionUnfolded = dimension;
	numChromosomes = dimension - 1;
	initStorage();
	fillFrequencies();
}

void SFS::initStorage(){
	clear();
	sfs = new float[dimension];
	sfsCumulative = new float[dimension];
	sfsFrequencies = new float[dimension];
	storageInitialized = true;
}

void SFS::clear(){
	if(storageInitialized){
		delete[] sfs;
		delete[] sfsCumulative;
		delete[] sfsFrequencies;
		storageInitialized = false;
	}
}

void SFS::fillFrequencies(){
	sfsFrequencies[0] = 0.0;
	for(int i=1; i<dimension; ++i){
		sfsFrequencies[i] = (float) i / (float) numChromosomes;
	}
}

void SFS::fillCumulative(){
	sfsCumulative[0] = sfs[0];
	for(int i=1; i<dimension; ++i){
		sfsCumulative[i] = sfsCumulative[i-1] + sfs[i];
	}
	sfsCumulative[numChromosomes] = 1.0;
}

void SFS::writeToFile(std::string & filename, bool writeLog){
	std::ofstream out(filename.c_str());
	if(!out)
		throw "Failed to open file '" + filename + "' for writing!";

	if(writeLog){
		out << log(sfs[0]);
		for(int i=1; i<dimension; ++i){
			out << " " << log(sfs[i]);
		}
	} else {
		out << sfs[0];
		for(int i=1; i<dimension; ++i){
			out << " " << sfs[i];
		}
	}
	out << "\n";
	out.close();
}

double SFS::calcLLOneSite(float* gl){
	double LL = 0.0;
	for(int i=0; i<dimension; ++i){
		LL += exp(gl[i]) * sfs[i];
	}

	return log(LL);
}

double SFS::getRandomFrequency(TRandomGenerator* randomGenerator){
	return sfsFrequencies[randomGenerator->pickOne(dimension, sfsCumulative)];
}

int SFS::getRandomAlleleCount(TRandomGenerator* randomGenerator){
	return randomGenerator->pickOne(dimension, sfsCumulative);
}

//--------------------------------------
//SFSfolded
//--------------------------------------
SFSfolded::SFSfolded(int numChr, float theta):SFS(numChr){
	storageInitialized = false;
	rescaled = false;
	initDimension(numChr  + 1);

	//generate sfs from theta
	float sum = 0;
	sfsFrequencies[0] = 0.0;
	for(int i=1; i<dimension; ++i){
		sfs[i] = theta / (float) i + theta / (float) (numChromosomes - i);
		sum += sfs[i];
		sfsFrequencies[i] = (float) i / (float) numChr;
	}
	if(sum > 1.0)
		throw "The choice of theta and sample size results in too many mutations in the SFS!";

	sfs[0] = 1.0 - sum;
	monoFrac = sfs[0];
	fillCumulative();
}

void SFSfolded::initDimension(int size){
	dimension = size;
	dimensionUnfolded = 2 * dimension - 1;
	numChromosomes = dimensionUnfolded - 1;
	initStorage();
	fillFrequencies();
}

double SFSfolded::calcLLOneSite(float* gl){
	double LL = 0.0;
	for(int i=0; i<dimension; ++i)
		LL += exp(gl[i]) * sfs[i];

	int j = dimension-1;
	for(int i=dimension; i<dimensionUnfolded; ++i, --j)
		LL += exp(gl[i]) * sfs[j];

	return log(LL);
}



