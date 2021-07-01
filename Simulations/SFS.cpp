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
SFS::SFS(const uint32_t & numChr){
	_rescaled = false;
	_initDimension(numChr + 1);
	monoFrac = 0.0;
};

SFS::SFS(const std::string & filename){
	_storageInitialized = false;
	std::ifstream file(filename.c_str());
	if(!file)
		throw "Failed to open SFS file '" + filename + "'!";
	std::vector<double> vec;
	coretools::str::fillContainerFromLineWhiteSpace(file, vec, true, true);

	//init dimension
	_initDimension(vec.size());

	//now store as fraction
	uint32_t i=0;
	float sum = 0;
	for(std::vector<double>::iterator it=vec.begin(); it!=vec.end(); ++it, ++i){
		sfs[i] = *it;
		sum += sfs[i];
	}
	for(i=0; i<dimension; ++i)
		sfs[i] = sfs[i] / sum;
	_rescaled = false;
	monoFrac = sfs[0];
	_fillCumulative();
};

SFS::SFS(SFS* other, const float & MonoFrac){
	//construct from other SFS, but rescale SFS to have a specific fraction of monomorphic sites
	_storageInitialized = false;
	_initDimension(other->dimension);

	//now copy and rescale
	_rescaled = true;
	monoFrac = MonoFrac;
	float sum = 0.0;
	sfs[0] = monoFrac;
	for(uint32_t i=1; i<dimension; ++i){
		sfs[i] = other->sfs[i];
		sum += sfs[i];
	}
	sum = sum / (1.0 - monoFrac);
	for(uint32_t i=1; i<dimension; ++i){
		sfs[i] = sfs[i] / sum;
	}
	_fillCumulative();
};

SFS::SFS(const uint32_t & numChr, const float & theta){
	_storageInitialized = false;
	_rescaled = false;
	_initDimension(numChr  + 1);

	//generate sfs from theta
	float sum = 0;
	sfsFrequencies[0] = 0.0;
	for(uint32_t i=1; i<dimension; ++i){
		sfs[i] = theta / (float) i;
		sum += sfs[i];
		sfsFrequencies[i] = (float) i / (float) numChr;
	}
	if(sum > 1.0)
		throw "The choice of theta and sample size results in too many mutations in the SFS!";

	sfs[0] = 1.0 - sum;
	monoFrac = sfs[0];
	_fillCumulative();
};

SFS::SFS(const uint32_t & numChr, const uint32_t & onlyThisBin){
	_storageInitialized = false;
	_rescaled = false;
	_initDimension(numChr  + 1);

	//set all to zero except this one bin
	for(uint32_t i=0; i<dimension; ++i){
		sfs[i] = 0.0;
	}
	sfs[onlyThisBin] = 1.0;
	monoFrac = sfs[0];
	_fillCumulative();
};

void SFS::_initDimension(const uint32_t & size){
	dimension = size;
	dimensionUnfolded = dimension;
	numChromosomes = dimension - 1;
	_initStorage();
	_fillFrequencies();
};

void SFS::_initStorage(){
	_clear();
	sfs = new float[dimension];
	sfsCumulative = new float[dimension];
	sfsFrequencies = new float[dimension];
	_storageInitialized = true;
};

void SFS::_clear(){
	if(_storageInitialized){
		delete[] sfs;
		delete[] sfsCumulative;
		delete[] sfsFrequencies;
		_storageInitialized = false;
	}
};

void SFS::_fillFrequencies(){
	sfsFrequencies[0] = 0.0;
	for(uint32_t i=1; i<dimension; ++i){
		sfsFrequencies[i] = (float) i / (float) numChromosomes;
	}
};

void SFS::_fillCumulative(){
	sfsCumulative[0] = sfs[0];
	for(uint32_t i=1; i<dimension; ++i){
		sfsCumulative[i] = sfsCumulative[i-1] + sfs[i];
	}
	sfsCumulative[numChromosomes] = 1.0;
};

void SFS::writeToFile(const std::string & filename, const bool & writeLog){
	std::ofstream out(filename.c_str());
	if(!out)
		throw "Failed to open file '" + filename + "' for writing!";

	if(writeLog){
		out << log(sfs[0]);
		for(uint32_t i=1; i<dimension; ++i){
			out << " " << log(sfs[i]);
		}
	} else {
		out << sfs[0];
		for(uint32_t i=1; i<dimension; ++i){
			out << " " << sfs[i];
		}
	}
	out << "\n";
	out.close();
};

double SFS::calcLLOneSite(float* gl){
	double LL = 0.0;
	for(uint32_t i=0; i<dimension; ++i){
		LL += exp(gl[i]) * sfs[i];
	}

	return log(LL);
};

double SFS::getRandomFrequency(coretools::TRandomGenerator* randomGenerator){
	return sfsFrequencies[randomGenerator->pickOne(dimension, sfsCumulative)];
};

uint32_t SFS::getRandomAlleleCount(coretools::TRandomGenerator* randomGenerator){
	return randomGenerator->pickOne(dimension, sfsCumulative);
};

//--------------------------------------
//SFSfolded
//--------------------------------------
SFSfolded::SFSfolded(const uint32_t & numChr, const float & theta):SFS(numChr){
	_storageInitialized = false;
	_rescaled = false;
	_initDimension(numChr  + 1);

	//generate sfs from theta
	float sum = 0;
	sfsFrequencies[0] = 0.0;
	for(uint32_t i=1; i<dimension; ++i){
		sfs[i] = theta / (float) i + theta / (float) (numChromosomes - i);
		sum += sfs[i];
		sfsFrequencies[i] = (float) i / (float) numChr;
	}
	if(sum > 1.0)
		throw "The choice of theta and sample size results in too many mutations in the SFS!";

	sfs[0] = 1.0 - sum;
	monoFrac = sfs[0];
	_fillCumulative();
};

void SFSfolded::_initDimension(const uint32_t & size){
	dimension = size;
	dimensionUnfolded = 2 * dimension - 1;
	numChromosomes = dimensionUnfolded - 1;
	_initStorage();
	_fillFrequencies();
};

double SFSfolded::calcLLOneSite(float* gl){
	double LL = 0.0;
	for(uint32_t i=0; i<dimension; ++i)
		LL += exp(gl[i]) * sfs[i];

	int j = dimension-1;
	for(uint32_t i=dimension; i<dimensionUnfolded; ++i, --j)
		LL += exp(gl[i]) * sfs[j];

	return log(LL);
};



