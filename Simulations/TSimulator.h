/*
 * TSimulator.h
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATOR_H_
#define TSIMULATOR_H_

#include <memory>
#include <string>
#include <vector>

#include "TSimulatorBamFiles.h"
#include "TSimulatorReference.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Types/probability.h"
#include "genometools/GenomePositions/TChromosomes.h"

#include "genometools/VCF/TVCFWriter.h"
#include "THaplotypeSimulator.h"
#include "TReadSimulators.h"

namespace Simulations {
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
	genometools::TChromosomes _chromosomes;
	std::unique_ptr<THaplotypeSimulator> _haploSimulator;

	virtual void _simulateAndWrite(const genometools::TChromosome &Chromosome, const TSimulatorHaplotypes &Haplotypes,
								   const TSimulatorReference &Reference, size_t avgDepth) = 0;

public:
	TSimulator(const std::string_view Method);
	void runSimulations();
	virtual ~TSimulator() = default;
};

//---------------------------------------------------------
// TBAMSimulator
//---------------------------------------------------------

class TBAMSimulator : public TSimulator {
protected:
	std::vector<TReadSimulators> _readSimulators; // one per sample
	TSimulatorBamFiles _bamFiles;

	void _simulateReadsForInd(const genometools::TChromosome &ThisChr, size_t Ind,
							  const std::vector<genometools::TwoBase> &Haplotypes, const TSimulatorReference &Reference,
							  TReadSimulators &ReadSimulator, size_t AvgDepth, BAM::TOutputBamFile &BamFile);

	void _simulateAndWrite(const genometools::TChromosome &Chromosome, const TSimulatorHaplotypes &Haplotypes,
						   const TSimulatorReference &Reference, size_t avgDepth) override;

public:
	TBAMSimulator(const std::string_view Method);
};

//---------------------------------------------------------
// TVCFSimulator
//---------------------------------------------------------

class TVCFSimulator : public TSimulator {
private:
	coretools::Probability _error{0.05};
	std::unique_ptr<genometools::TVCFWriter> _vcf;

protected:
	void _simulateAndWrite(const genometools::TChromosome &Chromosome, const TSimulatorHaplotypes &Haplotypes,
						   const TSimulatorReference &Reference, size_t avgDepth) override;

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
