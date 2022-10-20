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
#include "TSimulatorAuxiliaryTools.h"
#include "coretools/Math/TSubsamplePicker.h"
#include <functional>

namespace Simulations {
//--------------------------------
// Class to store and SFS
//--------------------------------

class SFS {
protected:
	SFS() = default;

	std::vector<size_t> _numChrPerPop;
	std::vector<size_t> _dimensions;
	size_t _numChr{};
	std::vector<double> sfs;
	std::vector<double> sfsCumulative;

	coretools::TSubsamplePicker _picker;

	void _setDerivedDiploid(const uint32_t l, TSimulatorHaplotypes & haplotypes, const size_t N, const size_t k, const size_t shift, const Base derived);
	void _setDerivedHaploid(const uint32_t l, TSimulatorHaplotypes & haplotypes, const size_t N, const size_t k, const size_t shift, const Base derived);

	size_t _simulateSite(const uint32_t l, TSimulatorHaplotypes & haplotypes, const Base ancestral, const Base derived, std::function<void(const uint32_t, TSimulatorHaplotypes &, const size_t, const size_t, const size_t, const Base)> func);

public:


	SFS(const std::string &filename);
	SFS(const SFS &other, float MonoFrac);
	SFS(uint32_t numChr, float theta);
	SFS(uint32_t numChr, uint32_t onlyThisBin);
	virtual ~SFS() = default;

	virtual uint32_t numChromosomes() const noexcept { return _numChr; };
	double monoFrac() const noexcept { return sfs.front(); };
	void writeToFile(const std::string &filename, const bool &writeLog = false);

	size_t simulateSiteDiploid(const uint32_t l, TSimulatorHaplotypes & haplotypes, const Base ancestral, const Base derived); //return true if site was polymorphic
	size_t simulateSiteHaploid(const uint32_t l, TSimulatorHaplotypes & haplotypes, const Base ancestral, const Base derived); //return true if site was polymorphic

	virtual double calcLLOneSite(const std::vector<float> &gl);
};

/*
class SFSfolded : public SFS {
public:
	virtual uint32_t numChromosomes() const noexcept override { return 2 * sfs.size() - 2; };
	SFSfolded(const std::string &filename) : SFS(filename){};
	SFSfolded(const SFSfolded &other, float MonoFrac) : SFS(other, MonoFrac){};
	SFSfolded(uint32_t numChr, float theta);
	double calcLLOneSite(const std::vector<float> &gl) override;
};
*/
} // namespace Simulations
#endif /* SFS_H_ */
