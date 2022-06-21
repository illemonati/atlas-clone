/*
 * SFS.h
 *
 *  Created on: Apr 10, 2017
 *      Author: wegmannd
 */

#ifndef SFS_H_
#define SFS_H_

#include <stdint.h>
#include <string>
#include <vector>

namespace Simulations {
//--------------------------------
// Class to store and SFS
//--------------------------------
class SFS {
protected:
	void _fillFrequencies();
	void _fillCumulative();
	SFS() = default;

public:
	std::vector<float> sfs;
	std::vector<float> sfsCumulative;
	std::vector<float> sfsFrequencies;

	SFS(const std::string &filename);
	SFS(const SFS &other, float MonoFrac);
	SFS(uint32_t numChr, float theta);
	SFS(uint32_t numChr, uint32_t onlyThisBin);
	virtual ~SFS() = default;

	virtual uint32_t numChromosomes() const noexcept { return sfs.size() - 1; };
	float monoFrac() const noexcept { return sfs.front(); };
	void writeToFile(const std::string &filename, const bool &writeLog = false);
	double getRandomFrequency();
	uint32_t getRandomAlleleCount();
	virtual double calcLLOneSite(const std::vector<float> &gl);
};

class SFSfolded : public SFS {
public:
	virtual uint32_t numChromosomes() const noexcept override { return 2 * sfs.size() - 2; };
	SFSfolded(const std::string &filename) : SFS(filename){};
	SFSfolded(const SFSfolded &other, float MonoFrac) : SFS(other, MonoFrac){};
	SFSfolded(uint32_t numChr, float theta) :SFS(numChr, theta){}
	double calcLLOneSite(const std::vector<float> &gl) override;
};

} // namespace Simulations
#endif /* SFS_H_ */
