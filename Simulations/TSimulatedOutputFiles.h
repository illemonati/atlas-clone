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
//#include "TAlignment.h"


namespace Simulations{

    class TSimulatedOutputFiles{

    private:
        std::vector<FASTQ::TFastqFile> _fastqFiles;        //need to be a pointer to single files //Simulations::TSimulatedOutputFile
        //std::vector<BAM::TOutputBamFile> _bamFiles;
        std::vector<std::string> filesIndex;                           //keep track of readGroups FASTQfiles

    public:
        TSimulatedOutputFiles();

        //explicit TSimulatedOutputFiles(std::vector<TSimulatedOutputFile *> *files);

        //initialization of vector of files

        //methods to operate on files
        void openFile(Simulations::TSimulatedOutputFile file);
        void addFile(std::string fileName);

        FASTQ::TFastqFile &operator[](size_t i);
    };

};      //namespace Simulations

#endif

