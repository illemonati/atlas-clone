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

	const auto ending = coretools::str::readAfterLast(_fileName, '.');										//get file suffix
	if(ending != NULL){UERROR("Impossible to create file from: ", _fileName, ", due to wrong filetype");}		//if suffix exists, print error
	else{ _fileName += ".fastq" ; }																			//if nothing, add .fastq
								//you sure it returns NULL in case nothing is found???
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

void TFastqFile::sortRead(BAM::TAlignment &alignment){
	uint16_t tmpReadGroupID = alignment.readGroupId();

	if(alignmentQueues.empty()) {
		static std::queue<BAM::TAlignment> newReadGroupIDAlignmentQueue;     //needs to be static in order to be reachable from other functions
		alignmentQueues.push(&newReadGroupIDAlignmentQueue);			
	}

	
}

bool TFastqFile::exists(uint16_t readGroupID){
	return false;
}

//------------------------------------------------
// Public methods
//------------------------------------------------


void TFastqFile::open(std::string_view fileName){
	tmpFastqFile.open(_fileName);		//open from Tfile
}

void TFastqFile::close(){ 
	tmpFastqFile.close();				//Close from TFile
 }

void TFastqFile::writeAlignment(BAM::TAlignment &alignment){
	sortRead(alignment);

	_writeAlignment(alignment);
}

void TFastqFile::setHeader(std::string header){  _header = header; }

void TFastqFile::setHeader(){ _header = genericIdentifiers; }

};		//namespace FASTQ

#endif