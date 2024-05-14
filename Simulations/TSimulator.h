/*
 * TSimulator.h
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATOR_H_
#define TSIMULATOR_H_

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "TSimulatorBamFiles.h"
#include "TSimulatorReference.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Types/probability.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "genometools/GenotypeTypes.h"

#include "genometools/VCF/TVcfWriter.h"
#include "THaplotypeSimulator.h"
#include "TReadSimulators.h"

namespace genometools { class PhredIntProbability; }

namespace Simulations {
class TSimulatorBamFiles;

// TODO: add cross-contamination between samples or RGs? That would be easier to model contamination that the way it is
// done now as it would allow for contaminated reads to have different characteristsics.

//---------------------------------------------------------
// TSimulator
//---------------------------------------------------------

class TSimulator {
protected:
	std::string _outname;
	std::vector<size_t> _seqDepth; //depth per chromosome
	bool _writeTrueGenotypes;
	bool _writeVariantInvariantBedFiles;
	TSimulatorReference _reference;
	genometools::TChromosomes _chromosomes;

	std::unique_ptr<THaplotypeSimulator> _haploSimulator;

	virtual void _simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes, size_t avgDepth) = 0;

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
	// bam files
	std::vector<TReadSimulators> _readSimulators; // one per sample
	std::unique_ptr<TSimulatorBamFiles> _bamFiles;

	// read simulator
	void _initializeReadSimulator();

	// functions to simulate
	void _simulateReadsFromHaplotypes(const genometools::TChromosome &thisChr,
									  const std::array<std::vector<genometools::Base>, 2>& haplotypes,
									  TReadSimulators &readSimulator, size_t avgDepth, BAM::TOutputBamFile &bamFile,
									  const std::string &extraProgressText);

	// simulate reads and write bam files
	void _simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes, size_t avgDepth) override;

public:
	TBAMSimulator(const std::string &method);
	~TBAMSimulator() { _bamFiles->close(); }
};

//---------------------------------------------------------
// TVCFSimulator
//---------------------------------------------------------

class TVCFSimulator : public TSimulator {
private:
	coretools::Probability _error{0.05};
	std::unique_ptr<genometools::TVcfWriter> _vcf;

protected:
	void _simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes, size_t avgDepth) override;

public:
	TVCFSimulator(const std::string &method);
};

struct TSimulationRunner {
	void run() {
		using coretools::instances::parameters;
		using coretools::instances::logfile;
		// initialize simulator
		auto method = parameters().get<std::string>("type", "one");

		if (parameters().exists("vcf")) {
			logfile().startIndent("Simulating VCF Files:");
			auto simulator = TVCFSimulator{method};
			simulator.runSimulations();
		} else { // default: BAM simulator
			logfile().startIndent("Simulating BAM Files:");
			auto simulator = TBAMSimulator{method};
			simulator.runSimulations();
		}

		// clean up
		logfile().endIndent();
	}
};

} // namespace Simulations

#endif /* TSIMULATOR_H_ */
