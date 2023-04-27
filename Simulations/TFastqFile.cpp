/*
 * TFastqFile.cpp
 *
 *  Created on: Marc 21, 2023
 *      Author: Michael Jopiti
 */


#ifndef TFASTQFILE_CPP
#define TFASTQFILE_CPP

#include <string>
#include <vector>
#include <stdlib.h>
#include "TFastqFile.h"
#include "genometools/PhredProbabilityTypes.h"
#include "TAlignment.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Math/counters.h"
#include "TAlignment.h"
#include "TReadSimulator.h"
#include "coretools/Main/TLog.h"
#include "coretools/Files/TWriter.h"
#include "coretools/Files/TReader.h"
#include "coretools/Files/TFile.h"
#include "coretools/Strings/stringFunctions.h" // for toString, readAfterLast, stringStartsWith
#include "TSimulatedOutputFile.h"

using coretools::instances::logfile; 	//used to write log file

namespace FASTQ{


//------------------------------------------------
// Methods to operate on Fastq files
//------------------------------------------------

    TFastqFile::TFastqFile(std::string_view fileName){
        _fileName = fileName;
        static coretools::TOutputFile txtFile(fileName);
    }


//------------------------------------------------
// Private methods
//------------------------------------------------

    void TFastqFile::_writeAlignment(const BAM::TAlignment &alignment){
        //alignment.flags
        //alignment.isEmpty
        //alignment.readGroupId
        //uint16_t readGroupId() const { return _readGroupID; };
        //alignment.refID
        //constexpr size_t refID() const noexcept { return _refID; };
        //alignment.qualities
        //std::vector<uint32_t> _qualities;
        //alignment.binQualityScoresIllumina
        //alignment.length
        //alignment.isPaired
        //alignment.mateRefID
        //alignment.clear
        //alignment.end


//takes alignment sequence and qualities and writes it in the file

    }

//------------------------------------------------
// Public methods
//------------------------------------------------


    void TFastqFile::open(std::string_view filename){
        //implement txtFile open()

    }

    void TFastqFile::close(){
        //implement txtFile close()
    }

    void TFastqFile::writeAlignment(const BAM::TAlignment &alignment){

        _writeAlignment(alignment);
    }

    void writeAlignmentLater(const BAM::TAlignment &alignment){
        //just to test if this method is the problem to lack of vtable
    }

    void TFastqFile::setHeader(std::string header){  _header = header; }

    void TFastqFile::setHeader(){}

    std::string_view TFastqFile::getName() {
        return _fileName;
    }

};		//namespace FASTQ

#endif