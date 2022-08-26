/*
 * THaplotypeSimulator.h
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#ifndef THAPLOTYPESIMULATOR_H_
#define THAPLOTYPESIMULATOR_H_

#include <stddef.h>
#include <stdint.h>
#include <array>
#include <memory>
#include <string>
#include <vector>

#include "GenotypeTypes.h"
#include "SFS.h"
#include "TFile.h"
#include "TGenotypeData.h"
#include "TSimulatorAuxiliaryTools.h"
#include "TStrongArray.h"

namespace genometools { class TChromosome; }
namespace genometools { class TChromosomes; }

namespace Simulations {

class THaplotypeSimulator {
protected:
	coretools::TStrongArray<double, genometools::Base> _cumulRef;
	GenotypeLikelihoods::TBaseProbabilities _baseFreq;
	 coretools::TStrongArray<double, genometools::Base>_cumulBaseFreq;
	coretools::Probability _referenceDivergence;
	THaplotypeSimulator();
public:
	virtual ~THaplotypeSimulator() = default;
	virtual void simulateDiploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference, const genometools::TChromosome &chromosome) = 0;
	virtual void simulateHaploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference, const genometools::TChromosome &chromosome) = 0;

	[[nodiscard]] virtual bool simulatesBiallelic() const noexcept = 0;
	[[nodiscard]] virtual size_t sampleSize() const noexcept = 0;
};

//---------------------------------------------------------
// TSimulatorOneIndividual
//---------------------------------------------------------
class TSimulatorOne : public THaplotypeSimulator {
private:
	std::vector<double> _thetas;

public:
	TSimulatorOne(size_t nChoromosomes);
	static inline const std::string name = "one";

	void simulateDiploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
					const genometools::TChromosome &chromosome) override;
	void simulateHaploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
					const genometools::TChromosome &chromosome) override;

	[[nodiscard]] bool simulatesBiallelic() const noexcept override { return true; };
	[[nodiscard]] size_t sampleSize() const noexcept override {return 1;}
};

//---------------------------------------------------------
// TSimulatorPairOfIndividuals
//---------------------------------------------------------
class TSimulatorPair : public THaplotypeSimulator {
private:
	static constexpr std::array<std::array<size_t, 4>, 4> _orderLookup = {
	    {{0, 1, 2, 3}, {0, 1, 3, 2}, {1, 0, 2, 3}, {1, 0, 3, 2}}};
	std::vector<double> _phis;
	std::array<double, 9> _cumulGenoCaseFrequencies;
	std::array<std::vector<double>, 9> _cumulGenoCombinationFreq;
	std::array<std::vector<std::array<genometools::Base, 4>>, 9> _genoTrans;

	void _fillTables();

public:
	TSimulatorPair();
	static inline const std::string name = "pair";
	void simulateDiploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
					const genometools::TChromosome &chromosome) override;
	void simulateHaploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
					const genometools::TChromosome &chromosome) override;

	[[nodiscard]] bool simulatesBiallelic() const noexcept override { return false; };
	[[nodiscard]] size_t sampleSize() const noexcept override {return 2;}
};

//---------------------------------------------------------
// TSimulatorSFS
//---------------------------------------------------------
class TSimulatorSFS : public THaplotypeSimulator {
private:
	int _sampleSize;
	std::vector<std::unique_ptr<SFS>> _sfs;
	TSimulatorMutationtable _mutTable;

	void _initializeSFS(const genometools::TChromosomes& chromosomes, const std::vector<double> &thetas, bool folded);
	void _initializeSFS(const genometools::TChromosomes& chromosomes, const std::vector<std::string> &sfsFileNames, bool folded);
public:
	TSimulatorSFS(const genometools::TChromosomes& chromosomes);
	static inline const std::string name = "SFS";
	void simulateHaploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
					const genometools::TChromosome &chromosome) override;
	void simulateDiploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
					const genometools::TChromosome &chromosome) override;

	[[nodiscard]] bool simulatesBiallelic() const noexcept override { return true; };
	[[nodiscard]] size_t sampleSize() const noexcept override {return _sampleSize;}
};

//---------------------------------------------------------
// TSimulatorHardyWeinberg
//---------------------------------------------------------
struct TSimulatorHWSite {
	bool isPolymorphic;
	genometools::Base reference;
	genometools::Base alternative;
	double f;
};

class TSimulatorHW : public THaplotypeSimulator {
private:
	int _sampleSize;
	double _fracPoly, _alpha, _beta, _F;
	double _cumulGenoProb[3];
	TSimulatorMutationtable _mutTable;
	bool _writeTrueAlleleFreq = false;
	coretools::TOutputFile _trueFreqFile;

	void _fillCumulGenoProb(double f);
	void _simulateSite(TSimulatorHWSite &site, TSimulatorReference &reference, const std::string &chr, uint64_t pos);
	void _fillhaplotypesMonomoprhic(TSimulatorHaplotypes &haplotypes, uint64_t locus, const TSimulatorHWSite &site);

public:
	TSimulatorHW();
	static inline const std::string name = "HW";
	void simulateHaploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
					const genometools::TChromosome &chromosome) override;
	void simulateDiploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
	                                const genometools::TChromosome &chromosome) override;

	[[nodiscard]] bool simulatesBiallelic() const noexcept override { return true; };
	[[nodiscard]] size_t sampleSize() const noexcept override {return _sampleSize;}
};
} // namespace Simulations

#endif
