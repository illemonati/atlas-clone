/*
 * SFS.h
 *
 *  Created on: Apr 10, 2017
 *      Author: wegmannd
 */

#ifndef SFS_H_
#define SFS_H_

#include <string>
#include <vector>

#include "TSimulatorHaplotypes.h"
#include "coretools/Main/TRandomPicker.h"
#include "coretools/Math/TSubsamplePicker.h"

namespace Simulations {
//--------------------------------
// Class to store and SFS
//--------------------------------

class SFS {
private:
	std::vector<size_t> _numChrPerPop;
	std::vector<size_t> _dimensions;
	size_t _numChr{};
	std::vector<double> _sfs;
	coretools::TRandomPicker _sfsPicker;

	coretools::TSubsamplePicker _picker;

	void _setDerivedDiploid(size_t l, TSimulatorHaplotypes & haplotypes, size_t N, size_t k, size_t shift, genometools::Base derived);
	void _setDerivedHaploid(size_t l, TSimulatorHaplotypes & haplotypes, size_t N, size_t k, size_t shift, genometools::Base derived);

	size_t _simulateSite(size_t l, TSimulatorHaplotypes & haplotypes, genometools::Base ancestral, genometools::Base derived, std::function<void(size_t, TSimulatorHaplotypes &, size_t, size_t, size_t, genometools::Base)> func);

public:
	SFS(const std::string &filename);
	SFS(const SFS &other, double MonoFrac);
	SFS(size_t numChr, double theta);
	SFS(size_t numChr, size_t onlyThisBin);

	size_t numChromosomes() const noexcept { return _numChr; };
	double monoFrac() const noexcept { return _sfs.front(); };
	void writeToFile(const std::string& filename, bool writeLog = false) const;

	size_t simulateSiteDiploid(size_t l, TSimulatorHaplotypes & haplotypes, genometools::Base ancestral, genometools::Base derived); //return true if site was polymorphic
	size_t simulateSiteHaploid(size_t l, TSimulatorHaplotypes & haplotypes, genometools::Base ancestral, genometools::Base derived); //return true if site was polymorphic

	double calcLLOneSite(const std::vector<double> &gl);
};
} // namespace Simulations
#endif /* SFS_H_ */
