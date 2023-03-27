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


using coretools::instances::logfile; //used to write log file

namespace FASTQ{

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

void TFastqFile::groupReads(){

}

//------------------------------------------------
// Public methods
//------------------------------------------------


void TFastqFile::open(std::string_view fileName){
    _fileName = fileName;

    logfile().list("Opening Fastq file '", _fileName, "'.");
	if (!_bamReader.Open(_fileName))                                    //is _bamReader unique per fastq file?
		UERROR("Failed to open BAM file '", _fileName, "'!");

	//load or create index file
	const std::string fnIndex1 = std::string(fileName).append(".fastqi");
	fileName.remove_suffix(4);
	const std::string fnIndex2 = std::string(fileName).append(".fastqi");
	if (std::filesystem::exists(fnIndex1)) {
		logfile().list("Opening Fastq index file '", fnIndex1, "'.");
		if(!_bamReader.OpenIndex(fnIndex1))
			UERROR("Failed to open Fastq index file '", fnIndex1, "'!");
	}
	else if (std::filesystem::exists(fnIndex2)) {
		logfile().list("Opening Fastq index file '", fnIndex2, "'.");
		if (!_bamReader.OpenIndex(fnIndex2))
			UERROR("Failed to open Fastq index file '", fnIndex2, "'!");
	} else {
		logfile().list("Creating Fastq index file '", fnIndex1, "'.");
		if (!_bamReader.CreateIndex())
			UERROR("Failed to create Fastq index file '", fnIndex1, "'!");
	}

    /**
     * write all the starting informations before writing reads and respective qualities
    */

    _open = true;
    //open();
}

void TFastqFile::close(){

    _open = false;
}

void TFastqFile::writeAlignment(BAM::TAlignment &alignment){

}

void TFastqFile::setHeader(std::string header){ 
    _header = header;
}

void TFastqFile::setHeader(){
    _header = genericIdentifiers;
}

}

#endif