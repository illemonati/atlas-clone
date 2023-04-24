//
// Created by Michael Jopiti on 19.04.23.
//

#ifndef SIMULATEDOUTPUTFILES_CPP
#define SIMULATEDOUTPUTFILES_CPP

#include <string>
#include "TSimulatedOutputFiles.h"
#include "../coretools/core/coretools/Files/TFile.h"

namespace Simulations{

    TSimulatedOutputFiles::TSimulatedOutputFiles(){}

    void TSimulatedOutputFiles::openFile(Simulations::TSimulatedOutputFile file) {

    }

    void TSimulatedOutputFiles::addFile(std::string fileName) {         //da migliorare il search
        //indexing new file with its file name
        filesIndex.push_back(fileName);

        //creation of new Fastq File
        FASTQ::TFastqFile newFile(fileName);
        _fastqFiles.push_back(newFile);
    }

    FASTQ::TFastqFile &TSimulatedOutputFiles::operator[](size_t i) {
        if (i >= _fastqFiles.size()) UERROR("FASTQ file ", i, " does not exist!");
        return _fastqFiles[i];
    }

};  //namespace Simulations

#endif
