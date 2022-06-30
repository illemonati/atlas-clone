/*
 * TSimulator.h
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATOR_H_
#define TSIMULATOR_H_

#include <stddef.h>
#include <stdint.h>
#include <array>
#include <exception>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "GenotypeTypes.h"
#include "TChromosomes.h"
#include "TGlfMultiReader.h"
#include "THaplotypeSimulator.h"
#include "TLog.h"
#include "TParameters.h"
#include "TPostMortemDamage.h"
#include "TReadGroups.h"
#include "TSequencingErrorModels.h"
#include "TSimulatorAuxiliaryTools.h"
#include "TSimulatorRead.h"
#include "TStrongArray.h"
#include "TTask.h"
#include "probability.h"

namespace genometools { class PhredIntProbability; }

namespace Simulations {

// TODO: add cross-contamination between samples or RGs? That would be easier to model contamination that the way it is
// done now as it would allow for contaminated reads to have different characteristsics.

//---------------------------------------------------------
// TSimulator
//---------------------------------------------------------

class TSimulator {
protected:
	std::string _outname;
	double _seqDepth;
	bool _writeTrueGenotypes;
	bool _writeVariantInvariantBedFiles;
	TSimulatorReference _reference;
	genometools::TChromosomes _chromosomes;

	std::unique_ptr<THaplotypeSimulator> _haploSimulator;

	virtual void _simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes) = 0;

public:
	TSimulator(const std::string &method);
	void runSimulations();
	virtual ~TSimulator() = default;
};

//---------------------------------------------------------
// TBAMSimulator
//---------------------------------------------------------

class TBAMSimulator : public TSimulator {
protected:
	double _averageReadLength = 0;
	double _maxReadLength     = 0;

	// bam files
	std::unique_ptr<TSimulatorBamFiles> _bamFiles;

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
	void _initializeSCDistribution(const std::string &ParameterName, const std::string &DefaultValue,
	                               const std::string &Name,
	                               std::function<void(TSimulatorSingleEndRead &, std::string, int distNumber)> function);
	void _initializeDistribution(const std::string &ParameterName, const std::string &DefaultValue,
	                             const std::string &Name,
	                             std::function<void(TSimulatorSingleEndRead &, std::string)> function);
	void _initializePMD(const std::string &ParameterName, const std::string &Name);
	void _initializeQualityTransformations(const std::string &ParameterName, const std::string &Name);
	void _initializeContamination(bool &perReadGroup, std::map<std::string, double> &contaminationMap);
	void _initializeReadSimulator();
	void _initializeReadGroupFrequencies();

	// functions to simulate
	void _simulateReadsFromHaplotypes(const genometools::TChromosome &thisChr,
	                                  std::array<std::vector<genometools::Base>, 2> haplotypes,
	                                  TSimulatorBamFile &bamFile, const std::string &extraProgressText);

	// simulate reads and write bam files
	void _simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes) override;

public:
	TBAMSimulator(const std::string &method);
};

//---------------------------------------------------------
// TVCFWriterSimulation
//---------------------------------------------------------

class TVCFWriterSimulation : public GLF::TGlfMultiReaderVcf {
public:
	TVCFWriterSimulation(const std::string &Filename, const std::string &Source,
	                     const std::vector<std::string> &SampleNames, bool UsePhredScaledLikelihoods);

	void writeSite(const std::string &ChrName, uint32_t Position, genometools::PhredIntProbability VariantQuality,
	               genometools::Base RefAllele, genometools::Base AltAllele, size_t TotalDepth, bool IsDiploid,
	               const GLF::TMultiGLFDataOneAllelicCombination &GenotypeLikelihoods,
	               const std::vector<size_t> &Depths);
};

//---------------------------------------------------------
// TVCFSimulator
//---------------------------------------------------------

class TVCFSimulator : public TSimulator {
private:
	coretools::Probability _error = 0.05;
	std::unique_ptr<TVCFWriterSimulation> _vcf;

protected:
	GLF::TMultiGLFDataSampleOneAllelicCombination _calculateGenotypeLikelihoods(size_t NumRef, size_t NumAlt,
	                                                                            bool IsDiploid);
	size_t _simulateNumReadsWithReferenceAllele(genometools::Base a, genometools::Base b, genometools::Base Ref,
	                                            size_t Depth, bool IsDiploid);
	std::pair<size_t, GLF::TMultiGLFDataSampleOneAllelicCombination>
	_simulateDepthAndGTL(genometools::Base a, genometools::Base b, genometools::Base Ref, bool IsDiploid);
	void _simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes) override;
	std::pair<genometools::Base, genometools::Base>
	_findMajorMinorAllele(coretools::TStrongArray<size_t, genometools::Base, 4> AlleleCounts,
	                      genometools::Base RefAllele);

public:
	TVCFSimulator(const std::string &method);
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
		std::unique_ptr<THaplotypeSimulator> haploSimulator;
		auto method = parameters().getParameterWithDefault<std::string>("type", "one");

		if (parameters().parameterExists("vcf")) {
			auto simulator = TVCFSimulator{method};
			simulator.runSimulations();
		} else { // default: BAM simulator
			auto simulator = TBAMSimulator{method};
			simulator.runSimulations();
		}

		// clean up
		logfile().endIndent();
	}
};

} // namespace Simulations

#endif /* TSIMULATOR_H_ */
