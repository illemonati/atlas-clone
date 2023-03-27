/*
 * TSimulatedOutputFile.h
 *
 *  Created on: Marc 22, 2023
 *      Author: Michael Jopiti
 */

#ifndef TSIMULATEDOUTPUTFILE_H
#define TSIMULATEDOUTPUTFILE_H

#include <string>
#include "../BAM/TBamFile.h"
#include "TFastqFile.h"
#include "../coretools/core/coretools/Files/TOutputFile.h"

namespace Simulations{

class TSimulatedOutputFile{
    private:
        coretools::TOutputFile simulatedOutputFile;

    public:
        virtual void open(std::string fileName) = 0;
        virtual void close() = 0;
        virtual void _writeAlignment(BAM::TAlignment &alignment) = 0;

};

/**
 * Useless this part because I just need to extend TSimulatedOutputFIle class in Bam and in FASTQ
*/

/*
class TSimulatedBamOutputFile : public TSimulatedOutputfile{
    // friend BAM::TOutputBamFile;      useless to call it friend because it has to be a wrapper not using eachothers functions

    public:
        TSimulatedBamOutputFile();  //still need to insert constructor istances
        void open() override;
        void close() override;
        void _writeAlignment() override;
};

class TSimulatedFastqOutputFile : public TSimulatedOutputfile{
    //friend FASTQ::TOutputFASTQFile;       useless in the same way for the TOutputBamFile

    public:
        TSimulatedFastqOutputFile();      //still need to insert constructor istances     
        void open();
        void close() override;
        void _writeAlignment() override;
};

*/

}

#endif