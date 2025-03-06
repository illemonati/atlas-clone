#ifndef TSIMULATORHAPLOTYPES_H_
#define TSIMULATORHAPLOTYPES_H_

#include "coretools/Files/TOutputFile.h"
#include "genometools/Genotypes/TwoBases.h"
#include <cstddef>
#include <vector>

namespace Simulations {
class TSimulatorReference;

class TSimulatorHaplotypes {
private:
	std::vector<std::vector<genometools::TwoBase>> _haplotypes;
	coretools::TOutputFile _trueGenoVCF; // write true genotypes to VCF
public:
	TSimulatorHaplotypes(size_t NumIndividuals) : _haplotypes(NumIndividuals) {}
	void setLength(size_t Length) noexcept {
		for (auto &h : _haplotypes) h.resize(Length);
	}
	size_t length() const { return size() ? _haplotypes.front().size() : 0; }
	void openTrueGenotypeVCF(std::string_view Filename);
	void writeTrueGenotypes(std::string_view ChrName, const TSimulatorReference &Ref);
	size_t size() const noexcept { return _haplotypes.size(); };
	const std::vector<genometools::TwoBase>& operator[](size_t Ind) const noexcept {return _haplotypes[Ind];}
	std::vector<genometools::TwoBase>& operator[](size_t Ind) noexcept {return _haplotypes[Ind];}
	bool isPolymoprhic(size_t Pos) const noexcept;
};
}
#endif
