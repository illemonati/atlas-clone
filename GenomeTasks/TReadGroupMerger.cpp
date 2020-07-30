/*
 * TReadGroupMerger.cpp
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#include "TReadGroupMerger.h"

namespace GenomeTasks{

TReadGroupMerger::TReadGroupMerger(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_basic(Parameters, Logfile, RandomGenerator){
	BAM::TReadGroups& readGroups = _bamFile.readGroupsMutable();

	//read read groups to be merged
	std::string filename = Parameters.getParameterString("readGroups");
	_logfile->startIndent("Reading read groups to be merged from file '" + filename + "':");
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
		fillVectorFromLineWhiteSpace(file, vec, true);
		if(!vec.empty()){
			if(vec.size() < 2) throw "Wrong number of entries on line " + toString(lineNum) + " in file '" + filename + "'!";

			//create new read group
			uint16_t newId = readGroups.add(vec[0]).id;
			_logfile->startIndent("The following read groups will be merged into '" + vec[0] + "':");

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
				_logfile->list(vec[i]);
			}
			_logfile->endIndent();
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
		_logfile->startIndent("The following read groups will be kept as is:");
		for(auto& s : unaffectedReadGroups){
			_logfile->list(s);
		}
		_logfile->endIndent();
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
	_bamFile.printEndWithSummary();

	//close bam writer
	outBam.close(_logfile);
};

}; // end namespace
