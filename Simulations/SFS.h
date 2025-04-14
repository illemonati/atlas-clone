/*
 * SFS.h
 *
 *  Created on: Apr 10, 2017
 *      Author: wegmannd
 */

#ifndef SFS_H_
#define SFS_H_

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

	size_t _simulateSite(size_t l, TSimulatorHaplotypes & haplotypes, genometools::Base ancestral, genometools::Base derived, bool Haplo);

public:
	SFS(std::string_view filename);
	SFS(size_t numChr, double theta);

	size_t numChromosomes() const noexcept { return _numChr; };
	double monoFrac() const noexcept { return _sfs.front(); };
	void writeToFile(std::string_view Filename, bool WriteLog = false) const;

	size_t simulateSiteDiploid(size_t L, TSimulatorHaplotypes & Haplotypes, genometools::Base Ancestral, genometools::Base Derived); 
	size_t simulateSiteHaploid(size_t L, TSimulatorHaplotypes & Haplotypes, genometools::Base Ancestral, genometools::Base Derived);
};
} // namespace Simulations
#endif /* SFS_H_ */
