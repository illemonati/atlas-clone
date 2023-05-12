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
    TFastqFile::TFastqFile(std::string_view fileName){      //: _file(fileName)
        _fileName = fileName;
        //_file.open(fileName);
    }

//------------------------------------------------
// Private methods
//------------------------------------------------


    std::string TFastqFile::_newName(int readGroupId) {
        if(readGroupId == BAM::TReadGroups::noReadGroupId){
            _outputName = _fileName + "_noReadGroup.fastq";
        } else _outputName = _fileName + "_" + std::to_string(readGroupId) + ".fastq";

        return _outputName;
    }

    void TFastqFile::_writeAlignment(coretools::TOutputFile& file, const BAM::TAlignment &alignment){
        //@readGroupID:refID:flags
        file.writeln("@" + std::to_string(alignment.readGroupId())
                      + ":" + std::to_string(alignment.refID())
                      + ":" + std::to_string(alignment.flags()));
        file.writeln(alignment.sequence());
        file.writeln("+");
        file.writeln(alignment.qualities());
    }

    bool TFastqFile::_idExists(uint16_t id){
        for (int i = 0; i < _readGroups.size(); ++i) {
            if (id == _readGroups[i]) {
                _exists = i;
                return true;
            }
        }
        return false;
    }

//------------------------------------------------
// Public methods
//------------------------------------------------

    void TFastqFile::open(std::string_view filename){
        _open = true;
        //_file.open(filename);
    }

    void TFastqFile::close() {
        _open = false;
        //_file.close();
    }

    void TFastqFile::writeAlignment(const BAM::TAlignment &alignment){

        if (_files.size() == 0){
            _readGroups.push_back(alignment.readGroupId());
            _files.emplace_back(new coretools::TOutputFile(_newName(alignment.readGroupId())));
            _writeAlignment(*_files[0], alignment);
        }else if (_idExists(alignment.readGroupId())){
            _writeAlignment(*_files[_exists], alignment);
        } else{
            _readGroups.push_back(alignment.readGroupId());
            _files.emplace_back(new coretools::TOutputFile(_newName(alignment.readGroupId())));
            _writeAlignment(*_files[_files.size() - 1], alignment);
        }

    }

    std::string_view TFastqFile::getName() {
        return _fileName;
    }

};		//namespace FASTQ

#endif