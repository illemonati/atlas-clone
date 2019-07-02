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
	bool rescaled;
	bool storageInitialized;
	virtual void initDimension(int size);
	void initStorage();
	void clear();
	void fillFrequencies();
	void fillCumulative();

public:
	float* sfs;
	float* sfsCumulative;
	float* sfsFrequencies;
	int numChromosomes;
	int dimension;
	int dimensionUnfolded;
	float monoFrac;

	SFS(int numChr);
	SFS(std::string filename);
	SFS(SFS* other, float MonoFrac);
	SFS(int numChr, float theta);
	SFS(int numChr, int onlyThisBin);
	virtual ~SFS(){ clear(); };

	void writeToFile(std::string & filename, bool writeLog=false);
	virtual double calcLLOneSite(float* gl);
	double getRandomFrequency(TRandomGenerator* randomGenerator);
	int getRandomAlleleCount(TRandomGenerator* randomGenerator);
};

class SFSfolded:public SFS{
public:
	SFSfolded(std::string filename):SFS(filename){};
	SFSfolded(SFSfolded* other, float MonoFrac):SFS(other, MonoFrac){};
	SFSfolded(int numChr, float theta);
	~SFSfolded(){};
	void initDimension(int size);
	double calcLLOneSite(float* gl);
};




#endif /* SFS_H_ */
