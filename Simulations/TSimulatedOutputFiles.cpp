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

#include <iostream>

namespace Simulations{

    TSimulatedOutputFiles::TSimulatedOutputFiles(uint32_t NumFiles, const std::string & Outname,
                                                 std::vector<TReadSimulators> & ReadSimulators,
                                                 const genometools::TChromosomes &Chromosomes) {

        coretools::instances::logfile().startIndent("Preparing FASTQ files for output:");

        _files.reserve(NumFiles);
        if (NumFiles < 1) DEVERROR("Can not open less than one FASTQ file!");
        if (ReadSimulators.size() > 1 && ReadSimulators.size() != NumFiles) {
            DEVERROR("Number of read simulators does not match number of files!");
        }
        if (Chromosomes.size() < 1) UERROR("Can not open a FASTQ file without specified chromosomes!");

        //addFile();          //create n fastqFiles
        for (uint32_t i = 0; i < NumFiles; i++) {
            addFile(NumFiles, Outname);
        }

    }

    Simulations::TSimulatedOutputFile &TSimulatedOutputFiles::operator[](size_t i) {
        if (i >= _files.size()) UERROR("FASTQ file ", i, " does not exist!");
        return *_files[i];
    }

    //what is the purpose?
    void TSimulatedOutputFiles::openFastqFile(FASTQ::TFastqFile fastqFile) {
        fastqFile.close();      //ensure the file is closed
        fastqFile.open(fastqFile.getName());
        _files.push_back(&fastqFile);
    }

    void TSimulatedOutputFiles::openBamFile(BAM::TOutputBamFile bamFile) {
        //BAM::open(bamFile.getOutputName(), &bamFile);
    }

    void TSimulatedOutputFiles::addFile(int NumFiles, const std::string & Outname) {     //creates 10 files
        if (NumFiles == 1) {
            _files.emplace_back(new FASTQ::TFastqFile(Outname));
        } else{
            _files.emplace_back(new FASTQ::TFastqFile(newName(Outname)));
        }
    }

    std::string TSimulatedOutputFiles::newName(const std::string & Outname) {
        _outputFileName = Outname + "_ind" + std::to_string(_files.size() + 1) ;
        return _outputFileName;
    }

};  //namespace Simulations

#endif
