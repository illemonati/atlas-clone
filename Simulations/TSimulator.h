/*
 * TSimulator.h
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATOR_H_
#define TSIMULATOR_H_

#include <algorithm>
#include <filesystem>
#include <functional>
#include <math.h>
#include <memory>
#include <numeric>

#include "SFS.h"
#include "TFile.h"
#include "TSimulatorAuxiliaryTools.h"
#include "TSimulatorRead.h"
#include "TTask.h"
#include "algorithmsAndVectors.h"
#include "progressTools.h"
#include "stringFunctions.h"
#include "THaplotypeSimulator.h"

namespace Simulations {

// TODO: add cross-contamination between samples or RGs? That would be easier to model contamination that the way it is
// done now as it would allow for contaminated reads to have different characteristsics.


//---------------------------------------------------------
// TSimulator
//---------------------------------------------------------

class TSimulator {
protected:
	std::string _outname;

	int _sampleSize             = 0;
	double _referenceDivergence = 0;
	double _seqDepth            = 0;

	genometools::TChromosomes _chromosomes;
	bool _writeTrueGenotypes            = false;
	bool _writeVariantInvariantBedFiles = false;
	TSimulatorReference _reference;

	std::unique_ptr<THaplotypeSimulator> _haploSimulator;
public:
	virtual void runSimulations();
	virtual ~TSimulator() = default;
};

class TVCFSimulator : public TSimulator {
	void runSimulations() override;
};


class TBAMSimulator : public TSimulator {
protected:
	double _averageReadLength   = 0;
	double _maxReadLength       = 0;

	// simulation tools
	BAM::TReadGroups _readGroups;
	GenotypeLikelihoods::TPostMortemDamage _PMD;
	// BAM::TReadGroupMap _readGroupMap; //needed by recal REALLYY??????
	GenotypeLikelihoods::SequencingError::TModels _recal;

	// read simulator
	std::vector<std::unique_ptr<TSimulatorSingleEndRead>> _readSimulators;
	std::vector<double> _simGroupFrequencies;
	std::vector<double> _cumulSimGroupFrequenies;


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

	// functions to simulate
	void _simulateReadsFromHaplotypes(const genometools::TChromosome &thisChr, std::array<std::vector<genometools::Base>,2> haplotypes, TSimulatorBamFile &bamFile,
					  std::string extraProgressText);
public:
	TBAMSimulator();
	// functions to set general parameters

	void runSimulations() override;
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
