/*
 * SFS.h
 *
 *  Created on: Apr 10, 2017
 *      Author: wegmannd
 */

#ifndef SFS_H_
#define SFS_H_

#include <iostream>
#include "stringFunctions.h"
#include <math.h>
#include "TRandomGenerator.h"

//--------------------------------
//Class to store and SFS
//--------------------------------
class SFS{
protected:
	bool _rescaled;
	bool _storageInitialized;
	virtual void _initDimension(const uint32_t & size);
	void _initStorage();
	void _clear();
	void _fillFrequencies();
	void _fillCumulative();

public:
	float* sfs;
	float* sfsCumulative;
	float* sfsFrequencies;
	uint32_t numChromosomes;
	uint32_t dimension;
	uint32_t dimensionUnfolded;
	float monoFrac;

	SFS(const uint32_t & numChr);
	SFS(const std::string & filename);
	SFS(SFS* other, const float & MonoFrac);
	SFS(const uint32_t & numChr, const float & theta);
	SFS(const uint32_t & numChr, const uint32_t & onlyThisBin);
	virtual ~SFS(){ _clear(); };

	void writeToFile(const std::string & filename, const bool & writeLog=false);
	virtual double calcLLOneSite(float* gl);
	double getRandomFrequency(coretools::TRandomGenerator* randomGenerator);
	uint32_t getRandomAlleleCount(coretools::TRandomGenerator* randomGenerator);
};

class SFSfolded:public SFS{
public:
	SFSfolded(const std::string & filename):SFS(filename){};
	SFSfolded(SFSfolded* other, const float & MonoFrac):SFS(other, MonoFrac){};
	SFSfolded(const uint32_t & numChr, const float & theta);
	~SFSfolded(){};
	void _initDimension(const uint32_t & size);
	double calcLLOneSite(float* gl);
};




#endif /* SFS_H_ */
