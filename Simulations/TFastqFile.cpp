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
#include "coretools/Strings/toString.h"
#include "coretools/Strings/stringFunctions.h" // for toString, readAfterLast, stringStartsWith
#include "TSimulatedOutputFile.h"

#include <iostream>

using coretools::instances::logfile; 	//used to write log file

namespace FASTQ{

    TFastqFile::TFastqFile(std::string_view fileName){
        _fileName = fileName;
        static coretools::TOutputFile txtFile(fileName);
        _ptrFile = &txtFile;
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
    //takes alignment sequence and qualities and writes it in the file
        _ptrFile->writeln("@" + std::to_string(alignment.readGroupId())
                            + ":" + std::to_string(alignment.refID()));
        _ptrFile->writeln(alignment.sequence());
        _ptrFile->writeln(alignment.qualities());
        _ptrFile->writeln("+");
    }

    std::string_view TFastqFile::getName() {
        return _fileName;
    }

};		//namespace FASTQ

#endif