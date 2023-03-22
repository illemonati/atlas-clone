/*
 * TSimulatedOutputFile.h
 *
 *  Created on: Marc 22, 2023
 *      Author: Michael Jopiti
 */

#include <string>
#include "../BAM/TBamFile.h"
#include "TFastqFile.h"

namespace Simulations{

class TSimulatedOutputfile{

    public:
        TSimulatedOutputfile();

        virtual void close() = 0;
        virtual void _writeAlignment() = 0;

};

class TSimulatedBamOutputFile : public TSimulatedOutputfile{
    friend BAM::TOutputBamFile;

    public:
        TSimulatedBamOutputFile();  //still need to insert constructor istances
        void close() override;
        void _writeAlignment() override;
};

class TSimulatedFASTQFile : public TSimulatedOutputfile{
    friend FASTQ::TOutputFASTQFile;

    public:
        TSimulatedFASTQFile();      //still need to insert constructor istances     
        void close() override;
        void _writeAlignment() override;
};

}