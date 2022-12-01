/*
 * TReadGroupMerger.cpp
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#include "TReadGroupMerger.h"

#include <stddef.h>
#include <algorithm>
#include <cstdint>
#include <exception>
#include <ios>
#include <set>
#include <utility>

#include "TBamFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "TReadGroups.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenomeTasks{
using coretools::instances::logfile;
using coretools::instances::parameters;

TReadGroupMerger::TReadGroupMerger():TGenome_basic(){
	BAM::TReadGroups& readGroups = _bamFile.readGroupsMutable();

	//read read groups to be merged
	std::string filename = parameters().getParameter<std::string>("readGroups");
	logfile().startIndent("Reading read groups to be merged from file '" + filename + "':");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open file '" + filename + "!";

	//create map oldId -> new Id. Fill with identity.
	readGroupMap.resize(readGroups.size());
	for(size_t i=0; i<readGroups.size(); ++i){
		readGroupMap[i] = i;
	}

	//parse file and construct new read groups in new header object
	int lineNum = 0;
	std::vector<std::string> vec;
	std::set<std::string> readGroupsMerged;
	while(file.good() && !file.eof()){
		++lineNum;
		coretools::str::fillContainerFromLineWhiteSpace(file, vec, true);
		if(!vec.empty()){
			if(vec.size() < 2) throw "Wrong number of entries on line " + coretools::str::toString(lineNum) + " in file '" + filename + "'!";

			//create new read group
			uint16_t newId = readGroups.add(vec[0]).id();
			logfile().startIndent("The following read groups will be merged into '" + vec[0] + "':");

			for(size_t i=1; i<vec.size(); ++i){
				//check for duplicates
				if(!readGroupsMerged.emplace(vec[i]).second){
					throw "Read group '" + vec[i] + "' is listed multiple times in file '" + filename + "'!";
				}

				uint16_t oldId = readGroups.getId(vec[i]);

				//set not to write to header
				readGroups.removeFromHeader(oldId);

				//update map
				readGroupMap[oldId] = newId;

				//report
				logfile().list(vec[i]);
			}
			logfile().endIndent();
		}
	}

	//report unaffected read groups
	std::vector<std::string> unaffectedReadGroups;
	for(size_t i=0; i<readGroups.size(); ++i){
		if(readGroupMap[i] ==  i){
			unaffectedReadGroups.emplace_back(readGroups.getName(i));
		}
	}

	if(unaffectedReadGroups.size() > 0){
		logfile().startIndent("The following read groups will be kept as is:");
		for(auto& s : unaffectedReadGroups){
			logfile().list(s);
		}
		logfile().endIndent();
	}
};

void TReadGroupMerger::mergeReadGroups(){
	//open a bam file for writing
	BAM::TOutputBamFile outBam;
	_openBamForWriting(_outputName + "_mergedRG.bam", outBam);

	//now parse through bam file and write alignments
	_bamFile.startProgressReporting();
	while(_bamFile.readNextAlignmentThatPassesFilters()){
		_bamFile.curSetNewReadGroup(readGroupMap[_bamFile.curReadGroupID()]);
		_bamFile.writeCurAlignment(outBam);

		//report
		_bamFile.printProgress();
	}
	_bamFile.printEndWithSummary(_outputName);

	//close bam writer
	outBam.close(&logfile());
};

}; // end namespace
