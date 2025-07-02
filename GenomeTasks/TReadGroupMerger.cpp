/*
 * TReadGroupMerger.cpp
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#include "TReadGroupMerger.h"
#include "TOutputBamFile.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TParameters.h"

namespace GenomeTasks{
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::user_assert;

TReadGroupMerger::TReadGroupMerger() {
	BAM::TReadGroups& readGroups = _genome.bamFile().readGroupsMutable();

	//read read groups to be merged
	const auto filename = parameters().get("readGroups");
	logfile().startIndent("Reading read groups to be merged from file '" + filename + "':");

	//create map oldId -> new Id. Fill with identity.
	_readGroupMap.resize(readGroups.size());
	for(size_t i=0; i<readGroups.size(); ++i){
		_readGroupMap[i] = i;
	}

	//parse file and construct new read groups in new header object
	std::set<std::string> readGroupsMerged;
	for (coretools::TInputFile file(filename.c_str(), coretools::FileType::NoHeader); !file.empty(); file.popFront()) {
		user_assert(file.numCols() >= 2, "Wrong number of entries on line ", file.curLine(), " in file '", filename, "'!");

		// create new read group
		uint16_t newId = readGroups.add(file.get(0)).id;
		logfile().startIndent("The following read groups will be merged into '", file.get(0), "':");

		for (size_t i = 1; i < file.numCols(); ++i) {
			// check for duplicates
			if (!readGroupsMerged.emplace(file.get(i)).second) {
				throw coretools::TUserError("Read group '", file.get(i), "' is listed multiple times in file '", filename, "'!");
			}

			const auto oldId = readGroups.getId(file.get(i));

			// set not to write to header
			readGroups.removeFromHeader(oldId);
			_readGroupMap[oldId] = newId;
			logfile().list(file.get(i));
		}
		logfile().endIndent();
	}

	//report unaffected read groups
	std::vector<std::string> unaffectedReadGroups;
	for(size_t i=0; i<readGroups.size(); ++i){
		if(_readGroupMap[i] ==  i){
			unaffectedReadGroups.push_back(readGroups.getName(i));
		}
	}

	if(unaffectedReadGroups.size() > 0){
		logfile().startIndent("The following read groups will be kept as is:");
		for(const auto& s : unaffectedReadGroups){
			logfile().list(s);
		}
		logfile().endIndent();
	}
}

void TReadGroupMerger::run(){
	//open a bam file for writing
	BAM::TOutputBamFile outBam(_genome.outputName() + "_mergedRG.bam", _genome.bamFile());

	//now parse through bam file and write alignments
	_genome.bamFile().startProgressReporting();
	while(_genome.bamFile().readNextAlignmentThatPassesFilters()){
		_genome.bamFile().curSetNewReadGroup(_readGroupMap[_genome.bamFile().curReadGroupID()]);
		_genome.bamFile().writeCurAlignment(outBam);

		_genome.bamFile().printProgress();
	}
	_genome.bamFile().printEndWithSummary(_genome.outputName());
}

} // end namespace
