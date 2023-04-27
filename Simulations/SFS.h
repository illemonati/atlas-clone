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
private:
	std::vector<size_t> _numChrPerPop;
	std::vector<size_t> _dimensions;
	size_t _numChr{};
	std::vector<double> sfs;
	std::vector<double> sfsCumulative;

	coretools::TSubsamplePicker _picker;

	void _setDerivedDiploid(size_t l, TSimulatorHaplotypes & haplotypes, size_t N, size_t k, size_t shift, Base derived);
	void _setDerivedHaploid(size_t l, TSimulatorHaplotypes & haplotypes, size_t N, size_t k, size_t shift, Base derived);

	size_t _simulateSite(size_t l, TSimulatorHaplotypes & haplotypes, Base ancestral, Base derived, std::function<void(size_t, TSimulatorHaplotypes &, size_t, size_t, size_t, Base)> func);

public:
	SFS(const std::string &filename);
	SFS(const SFS &other, double MonoFrac);
	SFS(size_t numChr, double theta);
	SFS(size_t numChr, size_t onlyThisBin);

	size_t numChromosomes() const noexcept { return _numChr; };
	double monoFrac() const noexcept { return sfs.front(); };
	void writeToFile(const std::string& filename, bool writeLog = false) const;

	size_t simulateSiteDiploid(size_t l, TSimulatorHaplotypes & haplotypes, Base ancestral, Base derived); //return true if site was polymorphic
	size_t simulateSiteHaploid(size_t l, TSimulatorHaplotypes & haplotypes, Base ancestral, Base derived); //return true if site was polymorphic

	double calcLLOneSite(const std::vector<double> &gl);
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
