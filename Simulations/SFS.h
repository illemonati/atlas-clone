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

namespace Simulations {
//--------------------------------
//Class to store and SFS
//--------------------------------
class SFS{
protected:
	void _fillFrequencies();
	void _fillCumulative();
	SFS() = default;
public:
	std::vector<float> sfs;
	std::vector<float> sfsCumulative;
	std::vector<float> sfsFrequencies;

	SFS(const std::string & filename);
	SFS(SFS* other, const float & MonoFrac);
	SFS(uint32_t numChr, const float & theta);
	SFS(uint32_t numChr, uint32_t onlyThisBin);
	virtual ~SFS() = default;

	virtual uint32_t numChromosomes() const noexcept {return sfs.size() - 1;};
	float monoFrac() const noexcept {return sfs.front();};
	void writeToFile(const std::string & filename, const bool & writeLog=false);
	virtual double calcLLOneSite(float* gl);
	double getRandomFrequency(coretools::TRandomGenerator* randomGenerator);
	uint32_t getRandomAlleleCount(coretools::TRandomGenerator* randomGenerator);
};

class SFSfolded:public SFS{
public:
	virtual uint32_t numChromosomes() const noexcept {return 2*sfs.size() - 2;};
	SFSfolded(const std::string & filename):SFS(filename){};
	SFSfolded(SFSfolded* other, const float & MonoFrac):SFS(other, MonoFrac){};
	SFSfolded(uint32_t numChr, const float & theta);
	~SFSfolded(){};
	double calcLLOneSite(float* gl);
};

} // namespace Simulations
#endif /* SFS_H_ */
