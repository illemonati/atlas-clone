/*
 * TReadGroupMerger.cpp
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#include "TReadGroupMerger.h"
#include "TOutputBamFile.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Files/TInputRcpp.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/splitters.h"

namespace GenomeTasks{
using coretools::instances::logfile;
using coretools::instances::parameters;

TReadGroupMerger::TReadGroupMerger() {
	BAM::TReadGroups& readGroups = _genome.bamFile().readGroupsMutable();

	//read read groups to be merged
	const auto filename = parameters().get("readGroups");
	logfile().startIndent("Reading read groups to be merged from file '" + filename + "':");

	//create map oldId -> new Id. Fill with identity.
	_readGroupMap.resize(readGroups.size());
	std::vector<bool> willChange(readGroups.size(), false);
	for(size_t i=0; i<readGroups.size(); ++i){
		_readGroupMap[i] = i;
	}

	//parse file and construct new read groups in new header object
	std::set<std::string> readGroupsMerged;
	coretools::TInputFile mergeFile(filename, coretools::FileType::Header);

	const auto idx = mergeFile.indices({"receiver", "donor"});
	for (; !mergeFile.empty(); mergeFile.popFront()) {
		// create new read group
		const auto newId  = readGroups.add(mergeFile.get(idx.front())).id;
		willChange[newId] = true;
		logfile().startIndent("The following read groups will be merged into '", mergeFile.get(0), "':");

		coretools::str::TSplitter spl(mergeFile.get(idx.back()), ',');

		for (auto s : spl) {
			// check for duplicates
			if (!readGroupsMerged.emplace(s).second) {
				throw coretools::TUserError("Read group '", s, "' is listed multiple times in file '", filename, "'!");
			}

			const auto oldId  = readGroups.getId(s);
			willChange[oldId] = true;

			// set not to write to header
			readGroups.removeFromHeader(oldId);
			_readGroupMap[oldId] = newId;
			logfile().list(s);
		}
		logfile().endIndent();
	}

	std::string s = "";
	for (size_t i = 0; i < readGroups.size(); ++i) {
		if (!willChange[i]) s += readGroups.getName(i) + ", ";
	}
	if (!s.empty()) {
		s.pop_back(); s.pop_back(); // remove last ", "
		logfile().list("The following read groups will be kept as is: ", s);
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
