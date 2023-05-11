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
#include "../BAM/TReadGroupInfo.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "TGlfMultiReader.h"
#include "THaplotypeSimulator.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "TPostMortemDamage.h"
#include "TReadGroups.h"
#include "TSimulatorAuxiliaryTools.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TTask.h"
#include "coretools/Types/probability.h"
#include "TReadSimulator.h"
#include "TReadSimulators.h"

#include "../BAM/TBamFile.h"
#include "TSimulatedOutputFile.h"
#include "TSimulatedOutputFiles.h"


namespace genometools { class PhredIntProbability; }
namespace Simulations { class TFastqBamSimulator; }

namespace Simulations {

// TODO: add cross-contamination between samples or RGs? That would be easier to model contamination that the way it is
// done now as it would allow for contaminated reads to have different characteristsics.

//---------------------------------------------------------
// TSimulator
//---------------------------------------------------------

    class TSimulator {
    protected:
        std::string _outname;
        std::vector<uint32_t> _seqDepth; //depth per chromosome
        bool _writeTrueGenotypes;
        bool _writeVariantInvariantBedFiles;
        TSimulatorReference _reference;
        genometools::TChromosomes _chromosomes;

        std::unique_ptr<THaplotypeSimulator> _haploSimulator;
        std::vector<TReadSimulators> _readSimulators; // one per sample


        virtual void _simulateAndWrite(const genometools::TChromosome &Chromosome,
                                       TSimulatorHaplotypes &Haplotypes, uint32_t avgDepth) = 0;

        //simulate reads to Fastq/bam file
        void _simulateReadsFromHaplotypes(const genometools::TChromosome &thisChr,
                                          std::array<std::vector<genometools::Base>, 2> haplotypes,
                                          Simulations::TReadSimulators &readSimulators, uint32_t avgDepth,
                                          Simulations::TSimulatedOutputFile &file, const std::string &extraProgressText);
        // read simulator
        void _initializeReadSimulator();

    public:
        TSimulator(const std::string &method);
        void runSimulations();
        virtual ~TSimulator() = default;
    };

//---------------------------------------------------------
// TBAMSimulator
//---------------------------------------------------------

    class TBAMSimulator : public TSimulator {

        friend class TFastqBamSimulator;

    protected:
        // bam files
        std::unique_ptr<TSimulatorBamFiles> _bamFiles;

        // simulate reads and write bam files
        void _simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes, uint32_t avgDepth) override;

    public:
        TBAMSimulator(const std::string &method);
        ~TBAMSimulator() = default;
    };


//-------------------------------------------
// FASTQ Simulator
//-------------------------------------------

    class TFastqSimulator : public TSimulator{

        friend class TFastqBamSimulator;

    private:
        //store created FastqFiles
        std::unique_ptr<Simulations::TSimulatedOutputFiles> _fastqFiles;

    public:
        //constructor / de constructor
        TFastqSimulator(const std::string &method);
        ~TFastqSimulator() = default;

        //overridden method from TSimulator (called by runSimulations() )
        void _simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes, uint32_t avgDepth) override;
    };

//-------------------------------------------
// FASTQ and BAM Simulator
//-------------------------------------------
class TFastqBamSimulator : public TSimulator{

    /**
     * create the fastq and bam simulators, but using the same seed.
     * they have to have the same simulated reads otherwise we have two different simulations
     * which are clearly not compatible and useful together
     * */

    private:
        //TFastqSimulator _fastqSimulator;
        //TBAMSimulator _bamSimulator;

        //store created FastqFiles
        std::unique_ptr<Simulations::TSimulatedOutputFiles> _fastqFiles;
        // bam files
        std::unique_ptr<TSimulatorBamFiles> _bamFiles;

        void _doubleSimulatedReadsFromHaplotypes(const genometools::TChromosome &thisChr,
                                                  std::array<std::vector<genometools::Base>, 2> haplotypes,
                                                  Simulations::TReadSimulators &readSimulators, uint32_t avgDepth,
                                                  Simulations::TSimulatedOutputFile &fastqFile,
                                                  BAM::TOutputBamFile &bamFile,
                                                  const std::string &extraProgressText);


    public:
        TFastqBamSimulator(const std::string &method);
        ~TFastqBamSimulator() = default;

        //overridden method from TSimulator (called by runSimulations() )
        void _simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes, uint32_t avgDepth) override;
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
        _simulateDepthAndGTL(genometools::Base a, genometools::Base b, genometools::Base Ref, bool IsDiploid, uint32_t avgDepth);
        void _simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes, uint32_t avgDepth) override;
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
            auto method = parameters().getParameterWithDefault<std::string>("type", "one");

            if (parameters().parameterExists("vcf")) {
                logfile().startIndent("Simulating VCF Files:");
                auto simulator = TVCFSimulator{method};
                simulator.runSimulations();
            } else if(parameters().parameterExists("fastq")){
                logfile().startIndent("Simulating FASTQ Files:");
                auto simulator = TFastqSimulator{method};
                simulator.runSimulations();
            } else if(parameters().parameterExists("fastq,bam") || parameters().parameterExists("bam,fastq")){
                logfile().startIndent("Simulating FASTQ and BAM files:");
                auto simulator = TFastqBamSimulator{method};
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
