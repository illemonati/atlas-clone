/*
 * TSimulator.h
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATOR_H_
#define TSIMULATOR_H_

#include "SFS.h"
#include "TFile.h"
#include "TSimulatorAuxiliaryTools.h"
#include "TSimulatorRead.h"
#include "TTask.h"
#include "algorithmsAndVectors.h"
#include "progressTools.h"
#include "stringFunctions.h"
#include <algorithm>
#include <filesystem>
#include <functional>
#include <math.h>
#include <memory>
#include <numeric>

namespace Simulations {

// TODO: add cross-contamination between samples or RGs? That would be easier to model contamination that the way it is
// done now as it would allow for contaminated reads to have different characteristsics.

using coretools::TLog;
using coretools::TParameters;
using coretools::TRandomGenerator;
using genometools::Base;

//---------------------------------------------------------
// TSimulator
//---------------------------------------------------------
class TSimulator {
protected:
	TLog *_logfile;
	TRandomGenerator *_randomGenerator;
	std::string _outname;

	// general simulation parameters
	int _sampleSize             = 0;
	double _referenceDivergence = 0;
	double _seqDepth            = 0;
	double _averageReadLength   = 0;
	double _maxReadLength       = 0;

	// chromosomes
	BAM::TChromosomes _chromosomes;
	bool _writeTrueGenotypes            = false;
	bool _writeVariantInvariantBedFiles = false;
	TSimulatorReference _referenceObj;

	// simulation tools
	BAM::TReadGroups _readGroups;
	GenotypeLikelihoods::TPostMortemDamage _PMD;
	// BAM::TReadGroupMap _readGroupMap; //needed by recal REALLYY??????
	GenotypeLikelihoods::TSequencingErrorModels _recal;

	// read simulator
	std::vector<std::unique_ptr<TSimulatorSingleEndRead>> _readSimulators;
	std::vector<double> _simGroupFrequencies;
	std::vector<double> _cumulSimGroupFrequenies;

	// helper tools
	std::array<double, 4> _cumulRef;
	GenotypeLikelihoods::TBaseProbabilities _baseFreq;
	std::array<double, 4> _cumulBaseFreq;
	bool _refInitialized = false;

	TSimulator(TLog *Logfile, TParameters &params, TRandomGenerator *RandomGenerator);

	// function to initialize read groups
	std::vector<std::string> _readSimInfoPerReadGroup(const std::string &Filename, const std::string &Column,
							  const std::string &Name);
	void _initializeReadGroup(const std::string &readLengthString, const BAM::TReadGroup &ReadGroup);
	void _initializeReadGroupsFromReadLengthDistribution(TParameters &params, const std::string &ParameterName,
							     const std::string &DefaultValue, const std::string &Name);
	void _initializeDistribution(TParameters &params, const std::string &ParameterName, const std::string &DefaultValue,
				     const std::string &Name,
				     std::function<void(TSimulatorSingleEndRead &, std::string)> function);
	void _initializePMD(TParameters &params, const std::string &ParameterName, const std::string &Name);
	void _initializeQualityTransformations(TParameters &params, const std::string &ParameterName,
					       const std::string &Name);
	void _initializeContamination(TParameters &params, bool &perReadGroup,
				      std::map<std::string, double> &contaminationMap);
	void _initializeChromosomes(TParameters &params);
	void _initializeReadSimulator(TParameters &params);
	void _initializeReadGroupFrequencies(TParameters &params);

	void _addToReadGroupVector(std::vector<std::string> &vec, const std::string &rg);
	void _addReadGroupsIfFile(const std::string &ParameterName, TParameters &Parameters, BAM::TReadGroups &ReadGroups);

	// functions to simulate
	Base _sampleBase(const std::array<double, 4> &cumulProbs);
	Base _mutateBase(const Base &base, const std::array<double, 4> &cumulProbs);
	virtual void _simulateHaplotypesDiploid(TSimulatorHaplotypes &, const BAM::TChromosome &) = 0;
	virtual void _simulateHaplotypesHaploid(TSimulatorHaplotypes &, const BAM::TChromosome &) = 0;
	void _simulateReadsFromHaplotypes(const BAM::TChromosome &thisChr, Base **haplotypes, TSimulatorBamFile &bamFile,
					  std::string extraProgressText);

public:
	virtual ~TSimulator() = default;
	// functions to set general parameters
	void setQualityDistribution(double mean, double sd, int maxQual);
	void setReadLength(std::string s);
	void setDepth(double depth) noexcept { _seqDepth = depth; };
	void setQualityTransformation(std::vector<double> &Betas);

	void runSimulations();
};

//---------------------------------------------------------
// TSimulatorOneIndividual
//---------------------------------------------------------
class TSimulatorOne : public TSimulator {
private:
	std::vector<double> thetas;

	void _simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;

public:
	TSimulatorOne(TLog *Logfile, TParameters &params, TRandomGenerator *RandomGenerator);
};

//---------------------------------------------------------
// TSimulatorPairOfIndividuals
//---------------------------------------------------------
class TSimulatorPair : public TSimulator {
private:
	static constexpr std::array<std::array<size_t, 4>, 4> orderLookup =
	{{{0, 1, 2, 3}, {0, 1, 3, 2}, {1, 0, 2, 3}, {1, 0, 3, 2}}};
	std::vector<double> phis;
	std::array<double, 9> cumulGenoCaseFrequencies;
	std::array<std::vector<double>, 9> cumulGenoCombinationFreq;
	std::array<std::vector<std::vector<Base>>, 9> genoTrans;

	void _fillTables();

	void _simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;

public:
	TSimulatorPair(TLog *Logfile, TParameters &params, TRandomGenerator *RandomGenerator);
};

//---------------------------------------------------------
// TSimulatorSFS
//---------------------------------------------------------
class TSimulatorSFS : public TSimulator {
private:
	std::vector<std::unique_ptr<SFS>> sfs;
	TSimulatorMutationtable mutTable;

	void _initializeSFS(const std::vector<double> &thetas);
	void _initializeSFS(const std::vector<std::string> &sfsFileNames, bool folded);
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;
	void _simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;

public:
	TSimulatorSFS(TLog *Logfile, TParameters &params, TRandomGenerator *RandomGenerator);
};

//---------------------------------------------------------
// TSimulatorHardyWeinberg
//---------------------------------------------------------
struct TSimulatorHW {
	bool isPolymorphic;
	genometools::Base reference;
	genometools::Base alternative;
	double f;
};

class TSimulatorHardyWeinberg : public TSimulator {
private:
	double fracPoly, alpha, beta, F;
	double cumulGenoProb[3];
	TSimulatorMutationtable mutTable;
	bool writeTrueAlleleFreq;
	coretools::TOutputFile trueFreqFile;

	void _fillCumulGenoProb(double f);
	void _simulateSite(TSimulatorHW &site, const std::string &chr, uint64_t pos);
	void _fillhaplotypesMonomoprhic(TSimulatorHaplotypes &haplotypes, uint64_t locus,
					TSimulatorHW &site);
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;
	void _simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;

public:
	TSimulatorHardyWeinberg(TLog *Logfile, TParameters &params, TRandomGenerator *RandomGenerator);
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_simulate : public coretools::TTask {
public:
	TTask_simulate() { _explanation = "Generating simulations"; };

	void run(TParameters &Parameters, TLog *Logfile) {
		// initialize simulator
		TSimulator *simulator;
		std::string method = Parameters.getParameterWithDefault<std::string>("type", "one");
		if (method == "one") {
			Logfile->startIndent("Simulating a single individual (parameter type=one):");
			simulator = new TSimulatorOne(Logfile, Parameters, _randomGenerator);
		} else if (method == "pair") {
			Logfile->startIndent("Simulating a pair of individual (parameter type=pair):");
			simulator = new TSimulatorPair(Logfile, Parameters, _randomGenerator);
		} else if (method == "SFS") {
			Logfile->startIndent("Simulating individuals from an SFS (parameter type=SFS):");
			simulator = new TSimulatorSFS(Logfile, Parameters, _randomGenerator);
		} else if (method == "HW") {
			Logfile->startIndent("Simulating a individuals under Hardy-Weinberg (parameter type=HW):");
			simulator = new TSimulatorHardyWeinberg(Logfile, Parameters, _randomGenerator);
		} else
			throw "Unknown simulation method '" + method + "'!";

		// now run simulations
		simulator->runSimulations();

		// clean up
		Logfile->endIndent();
		delete simulator;
	};
};

}; // namespace Simulations

#endif /* TSIMULATOR_H_ */
