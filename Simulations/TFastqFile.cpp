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


namespace FASTQ{

    using coretools::instances::logfile; //used to write log file

//------------------------------------------------
// Methods to operate on Fastq files
//------------------------------------------------

TFastqFile::TFastqFile(){
    _open = false;
}

//------------------------------------------------
// Private methods
//------------------------------------------------


void TFastqFile::_writeAlignment(BAM::TAlignment &alignment){

}

//------------------------------------------------
// Public methods
//------------------------------------------------


void TFastqFile::open(std::string fileName){
    _fileName = fileName;
    _open = true;
    //open();
}

/**
 * 
 * TBAMFile open method
 * 
 * void TOutputBamFile::open(const std::string Filename, const TBamFile & Original){
 *	open(Filename, Original.samHeader(), Original.chromosomes(), Original.readGroups());    
 *   };
*/

void TFastqFile::close(){

    _open = false;
}

void TFastqFile::writeAlignment(BAM::TAlignment &alignment){

}

void TFastqFile::setHeader(std::string header){ _header = header; }

}

#endif