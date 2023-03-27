/*
 * TBamFile.h
 *
 *  Created on: Marc 21, 2023
 *      Author: Michael Jopiti
 */

#ifndef TFASTQFILE_H_
#define TFASTQFILE_H_

#include <string>
#include <vector>
#include "genometools/PhredProbabilityTypes.h"
#include "TAlignment.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Math/counters.h"
#include "TAlignment.h"
#include "TReadSimulator.h"
#include "TSimulatedOutputFile.h"
#include "coretools/Main/TLog.h"

namespace FASTQ{

/*
 * 
 *  Create basic File set up for FASTQ file
 * 
*/

class TFastqFile : public Simulations::TSimulatedOutputFile{

    private:
    //istances regarding file handling
        std::string _fileName;
        std::string_view outputName;

    //generic inputs in case nothing is specified
        static constexpr std::string_view genericIdentifiers = "@FS001:001:0000000:1:1:0:0 1:Y:1:AAACCC";		
        static constexpr std::string_view genericFileName = "FASTQ_Simulation.fastq";

    //istances regarding file content (mainly for indexed file)
        int nbrReads;
        bool _open = false;
        int _maxSequenceLength = 80;
        std::string _header;
        std::string relativeQScoreValues;
        
    //istances regarding writing Fastq informations such as sequence and Phred quality score
        BAM::TAlignment *alignment;
        std::string _sequence;
        std::string _qualities;

    //methods operating on file content
        void _writeAlignment(BAM::TAlignment &alignment);

    public:
        TFastqFile();

        void open(std::string fileName) override;
        void close() override;
        void writeAlignment(BAM::TAlignment &alignment) override;

    //getters and setters
        void setHeader(std::string header);
        void setOutputName(std::string outputName);

        bool isOpen() const { return _open; };

};

}

#endif