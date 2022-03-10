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


//---------------------------------------------------------
// TSimulator
//---------------------------------------------------------
class TSimulator {
protected:
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
	TSimulatorReference _reference;

	// simulation tools
	BAM::TReadGroups _readGroups;
	GenotypeLikelihoods::TPostMortemDamage _PMD;
	// BAM::TReadGroupMap _readGroupMap; //needed by recal REALLYY??????
	GenotypeLikelihoods::SequencingError::TModels _recal;

	// read simulator
	std::vector<std::unique_ptr<TSimulatorSingleEndRead>> _readSimulators;
	std::vector<double> _simGroupFrequencies;
	std::vector<double> _cumulSimGroupFrequenies;

	// helper tools
	std::array<double, 4> _cumulRef;
	GenotypeLikelihoods::TBaseProbabilities _baseFreq;
	std::array<double, 4> _cumulBaseFreq;
	bool _refInitialized = false;

	TSimulator();

	// function to initialize read groups
	std::vector<std::string> _readSimInfoPerReadGroup(const std::string &Filename, const std::string &Column,
							  const std::string &Name);
	void _initializeReadGroup(const std::string &readLengthString, const BAM::TReadGroup &ReadGroup);
	void _initializeReadGroupsFromReadLengthDistribution(const std::string &ParameterName,
							     const std::string &DefaultValue, const std::string &Name);
	void _initializeDistribution(const std::string &ParameterName, const std::string &DefaultValue,
				     const std::string &Name,
				     std::function<void(TSimulatorSingleEndRead &, std::string)> function);
	void _initializePMD(const std::string &ParameterName, const std::string &Name);
	void _initializeQualityTransformations(const std::string &ParameterName,
					       const std::string &Name);
	void _initializeContamination(bool &perReadGroup,
				      std::map<std::string, double> &contaminationMap);
	void _initializeChromosomes();
	void _initializeReadSimulator();
	void _initializeReadGroupFrequencies();

	void _addToReadGroupVector(std::vector<std::string> &vec, const std::string &rg);
	void _addReadGroupsIfFile(const std::string &ParameterName, BAM::TReadGroups &ReadGroups);

	// functions to simulate
	genometools::Base _sampleBase(const std::array<double, 4> &cumulProbs);
	genometools::Base _mutateBase(genometools::Base base, const std::array<double, 4> &cumulProbs);
	virtual void _simulateHaplotypesDiploid(TSimulatorHaplotypes &, const BAM::TChromosome &) = 0;
	virtual void _simulateHaplotypesHaploid(TSimulatorHaplotypes &, const BAM::TChromosome &) = 0;
	void _simulateReadsFromHaplotypes(const BAM::TChromosome &thisChr, std::array<std::vector<genometools::Base>,2> haplotypes, TSimulatorBamFile &bamFile,
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
	std::vector<double> _thetas;

	void _simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;

public:
	TSimulatorOne();
};

//---------------------------------------------------------
// TSimulatorPairOfIndividuals
//---------------------------------------------------------
class TSimulatorPair : public TSimulator {
private:
	static constexpr std::array<std::array<size_t, 4>, 4> _orderLookup =
	{{{0, 1, 2, 3}, {0, 1, 3, 2}, {1, 0, 2, 3}, {1, 0, 3, 2}}};
	std::vector<double> _phis;
	std::array<double, 9> _cumulGenoCaseFrequencies;
	std::array<std::vector<double>, 9> _cumulGenoCombinationFreq;
	std::array<std::vector<std::array<genometools::Base,4>>, 9> _genoTrans;

	void _fillTables();
	void _simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;

public:
	TSimulatorPair();
};

//---------------------------------------------------------
// TSimulatorSFS
//---------------------------------------------------------
class TSimulatorSFS : public TSimulator {
private:
	std::vector<std::unique_ptr<SFS>> _sfs;
	TSimulatorMutationtable _mutTable;

	void _initializeSFS(const std::vector<double> &thetas);
	void _initializeSFS(const std::vector<std::string> &sfsFileNames, bool folded);
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;
	void _simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;

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

class TSimulatorHW : public TSimulator {
private:
	double _fracPoly, _alpha, _beta, _F;
	double _cumulGenoProb[3];
	TSimulatorMutationtable _mutTable;
	bool _writeTrueAlleleFreq = false;
	coretools::TOutputFile _trueFreqFile;

	void _fillCumulGenoProb(double f);
	void _simulateSite(TSimulatorHWSite &site, const std::string &chr, uint64_t pos);
	void _fillhaplotypesMonomoprhic(TSimulatorHaplotypes &haplotypes, uint64_t locus, const TSimulatorHWSite &site);
	void _simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;
	void _simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) override;

public:
	TSimulatorHW();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_simulate : public coretools::TTask {
public:
	TTask_simulate() { _explanation = "Generating simulations"; };

	void run() {
		using namespace coretools::instances;
		// initialize simulator
		std::unique_ptr<TSimulator> simulator;
		std::string method = parameters().getParameterWithDefault<std::string>("type", "one");
		if (method == "one") {
			logfile().startIndent("Simulating a single individual (parameter type=one):");
			simulator = std::make_unique<TSimulatorOne>();
		} else if (method == "pair") {
			logfile().startIndent("Simulating a pair of individual (parameter type=pair):");
			simulator = std::make_unique<TSimulatorPair>();
		} else if (method == "SFS") {
			logfile().startIndent("Simulating individuals from an SFS (parameter type=SFS):");
			simulator = std::make_unique<TSimulatorSFS>();
		} else if (method == "HW") {
			logfile().startIndent("Simulating a individuals under Hardy-Weinberg (parameter type=HW):");
			simulator = std::make_unique<TSimulatorHW>();
		} else
			throw "Unknown simulation method '" + method + "'!";

		// now run simulations
		simulator->runSimulations();

		// clean up
		logfile().endIndent();
	}
};

} // namespace Simulations

#endif /* TSIMULATOR_H_ */
