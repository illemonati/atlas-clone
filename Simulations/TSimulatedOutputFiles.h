/*
 * TSimulatedOutputFiles.h
 *
 *  Created on: Marc 22, 2023
 *      Author: Michael Jopiti
 */

/**
 *
 * Create structure for working with BAM or FASTQ files:
 *  - create vector populated from one of those two file types
 *  - create empty constructor
 *  - create two different open() methods (openFASTQ() and openBAM())
 *  - still have operators for the twos
 *
*/

#ifndef TSimulatedOutputFiles_H_
#define TSimulatedOutputFiles_H_

#include <vector>
#include "TSimulatedOutputFile.h"
#include "TFastqFile.h"
#include "TSimulatedOutputFile.h"
#include "../coretools/core/coretools/Main/TLog.h"
#include "TReadSimulators.h"

namespace Simulations{

    class TSimulatedOutputFiles{

    private:
        std::vector<Simulations::TSimulatedOutputFile*> _files;
        std::string _outputFileName;

        //generate new Name for the n individuals
        std::string newName(const std::string & Outname);

    public:
        TSimulatedOutputFiles(uint32_t NumFiles, const std::string & Outname, std::vector<TReadSimulators> & ReadSimulators,
                              const genometools::TChromosomes &Chromosomes);

        //methods to operate on files
        void openFastqFile(FASTQ::TFastqFile fastqFile);
        void openBamFile(BAM::TOutputBamFile bamFile);
        void addFile(int NumFile, const std::string & Outname);

        Simulations::TSimulatedOutputFile &operator[](size_t i);
    };

};      //namespace Simulations

#endif

