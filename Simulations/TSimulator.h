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
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Types/probability.h"
#include "genometools/GenomePositions/TChromosomes.h"

#include "genometools/TFastaReader.h"
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
	genometools::TFastaReader _fastaReader;
	std::vector<double> _seqDepth; //depth per chromosome
	bool _writeTrueGenotypes;
	bool _writeVariantInvariantBedFiles;
	genometools::TChromosomes _chromosomes;
	std::unique_ptr<THaplotypeSimulator> _haploSimulator;

	void _makeChromosomes();
	std::unique_ptr<THaplotypeSimulator> _makeHaploSimulator(std::string_view Method);
	virtual void _simulateAndWrite(const genometools::TChromosome &Chromosome, const TSimulatorHaplotypes &Haplotypes,
								   coretools::TView<genometools::Base> Reference, double avgDepth) = 0;

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
							  const std::vector<genometools::TwoBase> &Haplotypes, coretools::TView<genometools::Base> Reference,
							  TReadSimulators &ReadSimulator, double AvgDepth, BAM::TOutputBamFile &BamFile);

	void _simulateAndWrite(const genometools::TChromosome &Chromosome, const TSimulatorHaplotypes &Haplotypes,
						   coretools::TView<genometools::Base> Reference, double avgDepth) override;

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
						   coretools::TView<genometools::Base> Reference, double avgDepth) override;

public:
	TVCFSimulator(const std::string &method);
};

struct TSimulationRunner {
	void run() {
		using coretools::instances::parameters;
		using coretools::instances::logfile;
		// default type ist "one" if no fasta is given, HKY85 otherwise
		const auto type =
			parameters().get("type", parameters().exists("fasta") ? TSimulatorHKY85::name : TSimulatorOne::name);

		if (parameters().exists("vcf")) {
			logfile().startIndent("Simulating VCF Files:");
			auto simulator = TVCFSimulator{type};
			simulator.runSimulations();
		} else { // default: BAM simulator
			logfile().startIndent("Simulating BAM Files:");
			auto simulator = TBAMSimulator{type};
			simulator.runSimulations();
		}

		// clean up
		logfile().endIndent();
	}
};

} // namespace Simulations

#endif /* TSIMULATOR_H_ */
