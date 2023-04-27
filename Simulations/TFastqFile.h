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
#include <queue>
#include "genometools/PhredProbabilityTypes.h"
#include "TAlignment.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Math/counters.h"
#include "TAlignment.h"
#include "TReadSimulator.h"
#include "TSimulatedOutputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Files/TFile.h"
#include "coretools/Files/TOutputFile.h"

namespace Simulations { class TSimulatedOutputFile; }

namespace FASTQ{

    class TFastqFile : public Simulations::TSimulatedOutputFile {
        //create temporary file from public coretools::TFile, not good practice to do multiple inheritance

    private:
        //istances regarding file handling
        std::string _fileName;
        std::string_view outputName;
        //coretools::TOutputFile txtFile;

        //istances regarding file content (mainly for indexed file)
        bool _open = false;
        std::string _header;
        std::string _tmpRefID;
        uint16_t _readGroupID;

        //istances regarding writing Fastq information such as sequence, Phred quality score and queues for sorting alignments per group read
        BAM::TAlignment *alignment;
        std::string _sequence;
        std::string _qualities;

        //queue of pointers to queues of next-sorted alignments and related methods
        //static std::queue<std::queue<BAM::TAlignment>*> alignmentQueues;

        //we get alignment per alignment, no need for a queue?

        //methods operating on alignment content
        void _writeAlignment(const BAM::TAlignment &alignment);

    public:
        //constructor
        TFastqFile(std::string_view fileName);

        // File object
        //coretools::TOutputFile _file(fileName);

        //methods
        void open(std::string_view filename);
        void close();
        void writeAlignment(const BAM::TAlignment &alignment);
        void writeAlignmentLater(const BAM::TAlignment &alignment);

        //getters and setters
        void setHeader(std::string header);     //important to set Header because it is not inherited from TFile since it is not a table
        void setHeader();            //set with generic identifiers
        void setOutputName(std::string outputName);

        bool isOpen() const { return _open; };

        std::string_view getName();

    };

};      //namespace FASTQ

#endif