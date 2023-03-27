/*
 * TFastqFile.h
 *
 *  Created on: Marc 21, 2023
 *      Author: Michael Jopiti
 */

/*
 * 
 *  Create basic File set up for FASTQ file
 * 
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

class TFastqFile : public Simulations::TSimulatedOutputFile{

    private:
    //istances regarding file handling
        std::string _fileName;
        std::string_view outputName;
        BamTools::BamReader _bamReader;     //Added because I need to check if the file has been opened simultaneously

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
        void regroupSegments();

    public:
        TFastqFile();

        void open(std::string_view fileName) override;
        void close() override;
        void writeAlignment(BAM::TAlignment &alignment) override;

    //getters and setters
        void setHeader(std::string header);
        void setHeader();            //set with generic identifiers
        void setOutputName(std::string outputName);

        bool isOpen() const { return _open; };

};

}

#endif