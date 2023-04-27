//
// Created by Michael Jopiti on 19.04.23.
//
//
// Created by Michael Jopiti on 19.04.23.
//

#ifndef SIMULATEDOUTPUTFILES_CPP
#define SIMULATEDOUTPUTFILES_CPP

#include <string>
#include "TSimulatedOutputFiles.h"
#include "../coretools/core/coretools/Files/TFile.h"
#include "../coretools/core/coretools/Main/TLog.h"
#include "TReadSimulators.h"

namespace Simulations{

    TSimulatedOutputFiles::TSimulatedOutputFiles(uint32_t NumFiles, const std::string & Outname,
                                                 std::vector<TReadSimulators> & ReadSimulators,
                                                 const genometools::TChromosomes &Chromosomes){

        coretools::instances::logfile().startIndent("Preparing FASTQ files for output:");

        _fastqFiles.reserve(NumFiles);
        if (NumFiles < 1) DEVERROR("Can not open less than one FASTQ file!");
        if(ReadSimulators.size() > 1 && ReadSimulators.size() != NumFiles){
            DEVERROR("Number of read simulators does not match number of files!");
        }
        if (Chromosomes.size() < 1) UERROR("Can not open a FASTQ file without specified chromosomes!");

        addFile();          //create first fastqFile

    }

    void TSimulatedOutputFiles::openFile(Simulations::TSimulatedOutputFile file) {

    }

    void TSimulatedOutputFiles::addFile() {         //da migliorare il search

        //creation of new Fastq File
        FASTQ::TFastqFile newFile(newName());

        //indexing new file with its file name
        filesIndex.push_back(newFile.getName());

        //adding file to the store vector
        _fastqFiles.push_back(newFile);
    }

    FASTQ::TFastqFile &TSimulatedOutputFiles::operator[](size_t i) {
        if (i >= _fastqFiles.size()) UERROR("FASTQ file ", i, " does not exist!");
        return _fastqFiles[i];
    }

    std::string TSimulatedOutputFiles::newName() {
        _outputFileName = "FASTQ_Simulations_" + std::to_string(_fastqFiles.size()) + ".fastq";
        return _outputFileName;
    }

};  //namespace Simulations

#endif
