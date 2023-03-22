/*
 * TBamFile.h
 *
 *  Created on: Marc 21, 2023
 *      Author: Michael Jopiti
 */


/**
 * 
 *  Create basic File set up for FASTQ file
 * 
*/



// #ifndef SIMULATIONS_TFASTQFILE_H_
// #define SIMULATIONS_TFASTQFILE_H_

#include <string>
#include <vector>
#include "genometools/PhredProbabilityTypes.h"
#include "TAlignment.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Math/counters.h"
#include "TAlignment.h"
#include "TReadSimulator.h"

namespace FASTQ{

class TFastqFile{
    private:
        int nbrReads;
        bool _open{false};
        char _cPassFilter;
        int _maxSequenceLength = 80;
        std::string _header;
        std::string relativeQScoreValues;
        static constexpr std::string_view genericIdentifiers = "@FS001:001:0000000:1:1:0:0 1:Y:1:AAACCC";		//generic identifiers if not specified
        static constexpr std::string_view fileName = "FASTQ_Simulation_File.fastq";

    public:
        TFastqFile();   //do we want it?

};

class TOutputFASTQFile : public coretools::TOutputFile, public Simulations::TReadSimulator{     //shouldn't I be using some of TOutputFile?
    friend TFastqFile;

    private:
        std::string _fileName;
        bool _open;
    public:
        void _writeSequence();
    };

}