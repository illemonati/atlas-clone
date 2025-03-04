#ifndef TSIMULATORHAPLOTYPES_H_
#define TSIMULATORHAPLOTYPES_H_

#include "coretools/Files/gzstream.h"
#include "genometools/Genotypes/Base.h"
#include <cstddef>
#include <vector>
#include <array>

namespace Simulations {
class TSimulatorReference;

class TSimulatorHaplotypes {
private:
	size_t numInd;
	size_t _length = 0;
	std::vector<std::array<std::vector<genometools::Base>,2>> haplotypes;

	// write true genotypes to VCF
	gz::ogzstream trueGenoVCF;

	void allocateStorage();
public:
	TSimulatorHaplotypes(size_t NumIndividuals): numInd(NumIndividuals){};
	~TSimulatorHaplotypes() {
		if (trueGenoVCF) trueGenoVCF.close();
	}
	void setLength(size_t length) noexcept;
	size_t length() const { return _length; };
	void openTrueGenotypeVCF(std::string filename);
	const std::array<std::vector<genometools::Base>,2>& get(size_t i) const;
	void writeTrueGenotypes(const std::string &chrName, const TSimulatorReference &ref);
	size_t size() const noexcept { return numInd; };
	genometools::Base &operator()(size_t ind, size_t hap, size_t site) noexcept { return haplotypes[ind][hap][site]; };
	const genometools::Base &operator()(size_t ind, size_t hap, size_t site) const noexcept { return haplotypes[ind][hap][site]; };
	bool isPolymoprhic(size_t pos) const noexcept;
};
}
#endif
