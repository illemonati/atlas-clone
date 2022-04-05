/*
 * THaplotypeSimulator.h
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#ifndef THAPLOTYPESIMULATOR_H_
#define THAPLOTYPESIMULATOR_H_

#include <array>
#include <vector>
#include "TGenotypeData.h"
#include "TSimulatorAuxiliaryTools.h"

namespace Simulations {

class THaplotypeSimulator {
protected:
	std::array<double, 4> _cumulRef;
	GenotypeLikelihoods::TBaseProbabilities _baseFreq;
	std::array<double, 4> _cumulBaseFreq;
	int _sampleSize;
	double _referenceDivergence;
	THaplotypeSimulator();
public:
	virtual void simulateDiploid(TSimulatorHaplotypes &haplotypes, const genometools::TChromosome &chromosome);
	virtual void simulateHaploid(TSimulatorHaplotypes &haplotypes, const genometools::TChromosome &chromosome);
};

//---------------------------------------------------------
// TSimulatorOneIndividual
//---------------------------------------------------------
class TSimulatorOne : public THaplotypeSimulator {
private:
	std::vector<double> _thetas;

	void _simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes,
					const genometools::TChromosome &chromosome) override;
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes,
					const genometools::TChromosome &chromosome) override;

public:
	TSimulatorOne();
};

//---------------------------------------------------------
// TSimulatorPairOfIndividuals
//---------------------------------------------------------
class TSimulatorPair : public THaplotypeSimulator {
private:
	static constexpr std::array<std::array<size_t, 4>, 4> _orderLookup = {
	    {{0, 1, 2, 3}, {0, 1, 3, 2}, {1, 0, 2, 3}, {1, 0, 3, 2}}};
	std::vector<double> _phis std::array<double, 9> _cumulGenoCaseFrequencies;
	std::array<std::vector<double>, 9> _cumulGenoCombinationFreq;
	std::array<std::vector<std::array<genometools::Base, 4>>, 9> _genoTrans;

	void _fillTables();
	void _simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes,
					const genometools::TChromosome &chromosome) override;
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes,
					const genometools::TChromosome &chromosome) override;

public:
	TSimulatorPair();
};

//---------------------------------------------------------
// TSimulatorSFS
//---------------------------------------------------------
class TSimulatorSFS : public THaplotypeSimulator {
private:
	std::vector<std::unique_ptr<SFS>> _sfs;
	TSimulatorMutationtable _mutTable;

	void _initializeSFS(const std::vector<double> &thetas);
	void _initializeSFS(const std::vector<std::string> &sfsFileNames, bool folded);
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes,
					const genometools::TChromosome &chromosome) override;
	void _simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes,
					const genometools::TChromosome &chromosome) override;

public:
	TSimulatorSFS();
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
	double _fracPoly, _alpha, _beta, _F;
	double _cumulGenoProb[3];
	TSimulatorMutationtable _mutTable;
	bool _writeTrueAlleleFreq = false;
	coretools::TOutputFile _trueFreqFile;

	void _fillCumulGenoProb(double f);
	void _simulateSite(TSimulatorHWSite &site, const std::string &chr, uint64_t pos);
	void _fillhaplotypesMonomoprhic(TSimulatorHaplotypes &haplotypes, uint64_t locus, const TSimulatorHWSite &site);
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes,
					const genometools::TChromosome &chromosome) override;
	void _simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes,
	                                const genometools::TChromosome &chromosome) override;

public:
	TSimulatorHW();
};
} // namespace Simulations

#endif
