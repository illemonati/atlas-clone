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

//#include "TAlignment.h"


namespace Simulations{

    class TSimulatedOutputFiles{

    private:
        std::vector<FASTQ::TFastqFile> _fastqFiles;
        FASTQ::TFastqFile* _ptrNewFile;
        std::string _outputFileName;

        //std::vector<BAM::TOutputBamFile> _bamFiles;               WHAT TO DO WITH BAM?

        std::string newName(const std::string & Outname);

    public:
        TSimulatedOutputFiles(uint32_t NumFiles, const std::string & Outname, std::vector<TReadSimulators> & ReadSimulators,
                              const genometools::TChromosomes &Chromosomes);

        //explicit TSimulatedOutputFiles(std::vector<TSimulatedOutputFile *> *files);

        //initialization of vector of files

        //methods to operate on files
        void openFile(Simulations::TSimulatedOutputFile file);
        void addFile(int NumFile, const std::string & Outname);

        FASTQ::TFastqFile &operator[](size_t i);
    };

};      //namespace Simulations

#endif

