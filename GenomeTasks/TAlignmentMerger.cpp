/*
 * TAlignmentMerger.cpp
 *
 *  Created on: Oct 20, 2022
 *      Author: phaentu
 */

#include "TAlignmentMerger.h"

#include <math.h>
#include <stdlib.h>

#include <cstdint>
#include <exception>
#include <iostream>
#include <utility>
#include <vector>

#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "TAlignment.h"
#include "coretools/Files/TFile.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "TGenotypeData.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TReadGroups.h"
#include "TSequencedBase.h"
#include "coretools/Math/counters.h"
#include "coretools/Types/probability.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Types/strongTypes.h"

namespace GenomeTasks{

namespace AlignmentMerger{

using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using namespace GenomeTasks::BamFilter;
using namespace coretools::str;

//-----------------------------------------
// TAlignmentMergerReadGroupSettings
//-----------------------------------------
void TAlignmentMergerReadGroupSettings::initialize(BAM::TReadGroups & readGroups){
	//make sure set is empty
	_settings.clear();

	//check if we only merge without truncation
	if(parameters().parameterExists("pairedReadGroups")){
		std::string pairedRG = parameters().getParameter<std::string>("pairedReadGroups");
		if(pairedRG == "all"){
			//mark all as paired
			for(uint16_t rg=0; rg<readGroups.size(); ++rg){
				_settings.emplace(rg, paired, 0);
			}
		} else {
			//will merge a subset and treat others as unchanged
			std::vector<std::string> vec;
			fillContainerFromString(pairedRG, vec, ',');
			logfile().listFlush("Parsing read group names from parameter 'pairedReadGroups' ...");

			//get IDs
			std::set<uint16_t, std::less<>> pairedIds;
			for(auto n : vec){
				pairedIds.emplace(readGroups.getId(n));
			}

			//add as unchanged or paired
			for(uint16_t rg=0; rg<readGroups.size(); ++rg){
				if(pairedIds.find(rg) == pairedIds.end()){
					_settings.emplace(rg, unchanged, 0);
				} else {
					_settings.emplace(rg, paired, 0);
				}
			}
			logfile().done();
		}
		_printSummary();
	} else {
		//do we have to ignore read groups present in file?
		std::vector<std::string> vec;
		std::set<uint16_t> readGroupsToIgnore;
		if(parameters().parameterExists("ignoreReadGroups")){
			std::string ignoredReadGroupsFile = parameters().getParameter<std::string>("ignoreReadGroups");
			logfile().listFlush("Reading read groups to ignore from file '" + ignoredReadGroupsFile + "' ...");
			coretools::TInputFile in(ignoredReadGroupsFile, false);
			while(in.read(vec)){
				if(readGroups.readGroupExists(vec[0])){
					readGroupsToIgnore.insert(readGroups.getId(vec[0]));
				}
			}
			logfile().done();
		}

		//read file with read group settings
		std::string readGroupSettingsFile = parameters().getParameter<std::string>("readGroupSettings");
		logfile().listFlush("Reading read groups from file '" + readGroupSettingsFile + "' ...");
		coretools::TInputFile in(readGroupSettingsFile, {"readGroup", "seqType", "seqCycles"});

		//read settings file
		uint16_t numNotInUse = 0;
		while(in.read(vec)){
			//ignore "allReadGroups" from BAMD output
			if (vec[0] != "allReadGroups"){
				//get read group ID
				//read groups not in use will get a warning and be ignored
				uint16_t rgId = readGroups.getId(vec[0]);
				if(!readGroups.readGroupInUse(rgId)){
					++numNotInUse;
				} else {
					//parse max cycles
					uint16_t maxCycles = 0;
					if(vec[2] != "NA" && vec[2] != "-"){
						if(!stringContainsOnlyNumbers(vec[2])){
							throw "Error reading file '" + in.name() + "' on line " + toString(in.curLine()) + ": max cycles should be a number!";
						}
						maxCycles = convertString<int>(vec[2]);
					}

					//check for duplicate entries
					if(_settings.find(rgId) != _settings.end()){
						throw "Duplicate entry in file '" + in.name() + "' for read group '" + vec[0] + "'!";
					}

					//act based on seqeuncing type (second column). Ignored read groups will be marked as "unchanged"
					if(readGroupsToIgnore.find(rgId) != readGroupsToIgnore.end() || vec[1] == "unchanged"){
						_settings.emplace(rgId, unchanged, 0);
					} else if(vec[1] == "single"){
						if(maxCycles < 1){
							throw "Error reading file '" + in.name() + "' on line " + toString(in.curLine()) + ": max cycles must be > 0 for read groups of type 'single'!";
						}

						//add to settings and create truncated read group
						_settings.emplace(rgId, readGroups.addAlternativeRG(vec[0] + "_truncated", vec[0]).id(), single, maxCycles);

					} else if(vec[1] == "mixed"){
						if(maxCycles < 1){
							throw "Error reading file '" + in.name() + "' on line " + toString(in.curLine()) + ": max cycles must be > 0 for read groups of type 'mixed'!";
						}

						//add to settings and create truncated read group
						_settings.emplace(rgId, readGroups.addAlternativeRG(vec[0] + "_truncated", vec[0]).id(), mixed, maxCycles);

					} else if(vec[1] == "paired"){
						_settings.emplace(rgId, paired, 0);
					} else {
						throw "Error reading file '" + in.name() + "' on line " + toString(in.curLine()) + ": Unknown read group type '" + vec[1] + "'! Expected 'unchanged', 'single', 'mixed' or 'paired'.";
					}
				}
			}
		}

		//set missing read groups to "unchanged"
		for(uint16_t rg=0; rg < readGroups.size(); ++rg){
			if(readGroups.readGroupInUse(rg) && _settings.find(rg)==_settings.end()){
				_settings.emplace(rg, unchanged, 0);
			}
		}

		//summarize
		logfile().done();
		_printSummary();
		if(numNotInUse > 0){
			logfile().warning(std::to_string(numNotInUse) + " read group(s) are present in file '" + in.name() + "' but are marked as not in use!");
		}
	}
};

void TAlignmentMergerReadGroupSettings::_printSummary(){
	//count
	std::vector<uint16_t> counts(4, 0);
	for(auto& s : _settings){
		++counts[s.type];
	}

	//summarize
	if(counts[unchanged] > 0){ logfile().conclude(counts[unchanged], " read groups will remain unchanged."); }
	if(counts[single] > 0   ){ logfile().conclude(counts[single], " single-end read groups will be split."); }
	if(counts[mixed] > 0    ){ logfile().conclude(counts[mixed], " mixed read groups will be split and merged."); }
	if(counts[paired] > 0   ){ logfile().conclude(counts[paired], " paired read groups to be merged."); }
};

void TAlignmentMergerReadGroupSettings::setAllAsUnchanged(const BAM::TReadGroups & readGroups){
	_settings.clear();
	for(uint16_t rg=0; rg < readGroups.size(); ++rg){
		if(readGroups.readGroupInUse(rg)){
			_settings.emplace(rg, unchanged, 0);
		}
	}
};

bool TAlignmentMergerReadGroupSettings::needTruncation() const{
	for(auto& s : _settings){
		if(s.type == single || s.type == mixed)
			return true;
	}
	return false;
};

bool TAlignmentMergerReadGroupSettings::needsMerging() const{
	for(const auto& s : _settings){
		if(s.type == paired || s.type == mixed)
			return true;
	}
	return false;
};

ReadGroupType TAlignmentMergerReadGroupSettings::getType(const uint16_t readGroupId) const{
	return _settings.find(readGroupId)->type;
};

uint16_t TAlignmentMergerReadGroupSettings::getMaxCycles(const uint16_t readGroupId) const{
	return _settings.find(readGroupId)->maxCycles;
};

const TAlignmentMergerReadGroupSetting& TAlignmentMergerReadGroupSettings::getSettings(const uint16_t readGroupId) const{
	return *_settings.find(readGroupId);
};

//-----------------------------------------
// TAlignmentMergerType
//-----------------------------------------
uint16_t TAlignmentMerger::merge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	//check if reads overlap

	std::pair<uint32_t,bool> overlapLength = determineOverlapLength(alignment, mate);
	if (overlapLength.first > 0){
		//if the second read is being merged, it's position in the BAM-file as well as the position of the mate of the first read need to be adjusted to account for the length of the added softclips on left side
		if(!overlapLength.second){
			alignment.moveOnRef(alignment.position() + overlapLength.first);
			mate.setMateGenomicPosition(alignment);
		}
		alignment.merge(overlapLength.first, overlapLength.second);
	}
	return overlapLength.first;
};

std::pair<uint32_t,bool> TAlignmentMerger::determineOverlapLength(const BAM::TAlignment & alignment, const BAM::TAlignment & mate){
	if (alignment.cigar().lengthAligned() == 0 || mate.cigar().lengthAligned() == 0)
		return std::make_pair(0,true);
	if (alignment > mate){
		if (alignment.position()  < mate.position() + mate.cigar().lengthAligned()){
			uint32_t overlapLength = mate.position() + mate.cigar().lengthAligned() - alignment.position();
			bool alignmentFirst = false;
			return std::make_pair(overlapLength,alignmentFirst);
		} return std::make_pair(0,false);
	} else {
		if (mate.position()  < alignment.position() + alignment.cigar().lengthAligned()){
			uint32_t overlapLength = alignment.position() + alignment.cigar().lengthAligned() - mate.position();
			bool alignmentFirst = true;
			return std::make_pair(overlapLength,alignmentFirst);
		} return std::make_pair(0,true);
	}
}

// TAlignmentMergerType_randomRead
//---------------------------------
TAlignmentMerger_randomRead::TAlignmentMerger_randomRead():TAlignmentMerger(){};

uint16_t TAlignmentMerger_randomRead::merge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	if (randomGenerator().pickOneOfTwo())
		return TAlignmentMerger::merge(mate, alignment);
	else
		return TAlignmentMerger::merge(alignment, mate);
};

// TAlignmentMergerType_firstMate
//---------------------------------
TAlignmentMerger_firstMate::TAlignmentMerger_firstMate():TAlignmentMerger(){};

uint16_t TAlignmentMerger_firstMate::merge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	if (alignment.isSecondMate())
		return TAlignmentMerger::merge(mate, alignment);
	else
		return TAlignmentMerger::merge(alignment, mate);
};

// TAlignmentMergerType_secondMate
//---------------------------------
TAlignmentMerger_secondMate::TAlignmentMerger_secondMate():TAlignmentMerger(){};

uint16_t TAlignmentMerger_secondMate::merge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	if (alignment.isSecondMate())
		return TAlignmentMerger::merge(alignment, mate);
	else
		return TAlignmentMerger::merge(mate, alignment);
};

// TAlignmentMergerType_highestQuality
//---------------------------------
TAlignmentMerger_highestQuality::TAlignmentMerger_highestQuality():TAlignmentMerger(){};

uint16_t TAlignmentMerger_highestQuality::merge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	std::pair<uint32_t,bool> overlapLength = TAlignmentMerger::determineOverlapLength(alignment, mate);
	if (overlapLength.first != 0) {
		std::pair<genometools::PhredIntProbability,genometools::PhredIntProbability> minQuals = getMinQuals(alignment, mate, overlapLength);
		if (minQuals.first > minQuals.second){
			if(overlapLength.second){
				alignment.moveOnRef(alignment.position() + overlapLength.first);
				mate.setMateGenomicPosition(alignment);
			}
			mate.merge(overlapLength.first, overlapLength.second);
			return overlapLength.first;
		} else if (minQuals.second > minQuals.first){
			if(!overlapLength.second){
				alignment.moveOnRef(alignment.position() + overlapLength.first);
				mate.setMateGenomicPosition(alignment);
			}
			alignment.merge(overlapLength.first, overlapLength.second);
			return overlapLength.first;
		} else {
				if (randomGenerator().pickOneOfTwo()) {
					if(overlapLength.second){
						alignment.moveOnRef(alignment.position() + overlapLength.first);
						mate.setMateGenomicPosition(alignment);
					}
					mate.merge(overlapLength.first, overlapLength.second);
					return overlapLength.first;
				} else {
					if(!overlapLength.second){
						alignment.moveOnRef(alignment.position() + overlapLength.first);
						mate.setMateGenomicPosition(alignment);
					}
					alignment.merge(overlapLength.first, overlapLength.second);
					return overlapLength.first;
				}
		}
	} else
		return 0;
}

std::pair<genometools::PhredIntProbability,genometools::PhredIntProbability> TAlignmentMerger_highestQuality::getMinQuals(BAM::TAlignment & alignment, BAM::TAlignment & mate, std::pair<uint16_t,bool> overlapLength) const {
	if (overlapLength.second == true)
		return minQual(alignment,mate);
	else {
		std::pair<genometools::PhredIntProbability,genometools::PhredIntProbability> flippedResult = minQual(mate,alignment);
		return std::make_pair(flippedResult.second,flippedResult.first);
	}
}

std::pair<genometools::PhredIntProbability,genometools::PhredIntProbability> TAlignmentMerger_highestQuality::minQual(BAM::TAlignment & firstRead, BAM::TAlignment & secondRead) const {
	auto baseIterator = secondRead.begin();
	uint16_t internalPos = 0;
	genometools::PhredIntProbability secondReadMinQual = baseIterator->recalibratedQualityAsPhredInt;
	while(secondRead.positionInRef(internalPos).position() != firstRead.lastAlignedPositionWithRespectToRef().position() && secondRead.positionInRef(internalPos).position()  != secondRead.lastAlignedPositionWithRespectToRef().position()){
		if(secondRead.isAlignedAtInternalPos(internalPos)){
			if (baseIterator->recalibratedQualityAsPhredInt < secondReadMinQual)
				secondReadMinQual = baseIterator->recalibratedQualityAsPhredInt;
		}
		baseIterator++;
		internalPos++;
	}

	baseIterator = --firstRead.end();
	internalPos = firstRead.getLastInternalPos();
	genometools::PhredIntProbability firstReadMinQual = baseIterator->recalibratedQualityAsPhredInt;
	while (firstRead.positionInRef(internalPos).position() != secondRead.position() && firstRead.positionInRef(internalPos).position()  != firstRead.lastAlignedPositionWithRespectToRef().position()) {
		if(firstRead.isAlignedAtInternalPos(internalPos)){
			if (baseIterator->recalibratedQualityAsPhredInt < firstReadMinQual)
				firstReadMinQual = baseIterator->recalibratedQualityAsPhredInt;
		}
		baseIterator--;
		internalPos++;
	}
	return std::make_pair(firstReadMinQual,secondReadMinQual);
}


//-----------------------------------------
// TAlignmentSplitMerger
//-----------------------------------------
TAlignmentSplitMerger::TAlignmentSplitMerger() : TGenomeParsedWithAlignmentStorage() {
	//parse read group settings
	_rgSettings.initialize(_bamFile.readGroupsMutable());

	//allow for reads to exceed max cycle length?
	if(_rgSettings.needTruncation() || parameters().parameterExists("allowForLarger")){
		_allowForLarger = true;
		logfile().list("Adding single end reads that are longer than maxCycles to 'truncated' read group without throwing an error. (parameter 'allowForLarger')");
	} else {
		_allowForLarger = false;
	}

	//initialize merger, if needed
	if(_rgSettings.needsMerging()){
		_initializeMerger();
	}
};

void TAlignmentSplitMerger::_initializeMerger() {
	// check if keepAllReads is turned on
	// TODO: what is the basic set of filters needed?
	if(!_bamFile.improperPairsFilterEnabled()){
		logfile().warning("Improper pairs are kept but will not be merged!");
	}
	/*
	if(alignmentParser.getKeepAll()){
		logfile().warning("Undefined behavior when merging reads that do not pass default filters. Consider removing 'keepAllReads'");
	}


	//decide if we update quality score
	bool adaptQuality;
	if(parameters().parameterExists("updateQuality")){
		adaptQuality = true;
		logfile().list("Will update quality scores of preferred bases to reflect information from overlapping bases.");
	} else {
		adaptQuality = false;
		logfile().list("Will keep original quality scores of the preferred bases (use updateQuality to update quality scores).");
	}*/

	//set merging method
	//TODO: update wiki to reflect change in names
	std::string method = parameters().getParameterWithDefault<std::string>("mergingMethod", "randomRead");
	if(method == "none"){
		_merger = std::make_unique<TAlignmentMerger>();
		logfile().list("Merging method: no merging. (parameter 'mergingMethod')");
	} else if (method == "randomRead"){
		_merger = std::make_unique<TAlignmentMerger_randomRead>();
		logfile().list("Merging method: will keep random read for all overlapping positions. (parameter 'mergingMethod')");
	} else if(method == "highestQuality"){
		_merger = std::make_unique<TAlignmentMerger_highestQuality>();
		logfile().list("Merging method: will keep read with highest quality at overlapping positions. (parameter 'mergingMethod')");
	} else if(method == "firstMate"){
		_merger = std::make_unique<TAlignmentMerger_firstMate>();
		logfile().list("Merging method: will keep first mate at overlapping positions. (parameter 'mergingMethod')");
	} else if(method == "secondMate"){
		_merger = std::make_unique<TAlignmentMerger_secondMate>();
		logfile().list("Merging method: will keep second mate at overlapping positions. (parameter 'mergingMethod')");
	} else {
		throw "Unknown merging method " + method + "! Use 'none', 'firstMate', 'secondMate', 'randomRead' or 'highestQuality'.";
	}
};

void TAlignmentSplitMerger::_openBamFileForWriting(){
	TGenome_basic::_openBamForWriting(_outputName + "_splitMerged.bam", _outBam);
};

void TAlignmentSplitMerger::_handleMates(BAM::TAlignment & alignment, TAlignmentStorageSortedIterator mate){
	ReadGroupType type = _rgSettings.getType(alignment.readGroupId());

	if(type == single){
		throw "Paired reads found in single-end read group '" + _bamFile.readGroups().getName(alignment.readGroupId()) + "'! Is this a 'mixed' read group?";
	} else if(!alignment.isProperPair()){
		//not a proper pair: mark mate as as improper too
		mate->setAsNonProperPair();
		mate->makeReady();
	} else if(type == paired || type == mixed){
		//attempt merging: make sure alignments are parsed
		//Note: if we recalibrate, they were already parsed
		if(!alignment.isParsed()){
			alignment.parse();
		}

		//since mate position couldbe affected by merge: extract from storage and put back in
		BAM::TAlignment* mateAlignment = mate->stealAlignment();
		_alignmentStorage.erase(mate);		
		
		if(!mateAlignment->isParsed()){
			mateAlignment->parse();
		}

		_merger->merge(alignment, *mateAlignment);		
		addToContainer(_alignmentStorage, mateAlignment, true);
		}

	//add alignment to container
	addToContainer(_alignmentStorage, &alignment, true);
};

void TAlignmentSplitMerger::_handleSingle(BAM::TAlignment & alignment){
	const TAlignmentMergerReadGroupSetting& settings = _rgSettings.getSettings(alignment.readGroupId());

	if(settings.type == unchanged){
		//add as ready for writing
		addToContainer(_alignmentStorage, &alignment, true);
	} else if(settings.type == single || settings.type == mixed){
		//truncate
		if(!_allowForLarger && alignment.length() > settings.maxCycles){
			throw("Length of read " + alignment.name() + " is > max cycles for its read group (" + toString(settings.maxCycles) + ")! Use parameter 'allowForLarger' to ignore and put read in truncated read group.");
		} else if(alignment.length() >= settings.maxCycles){
			//add to truncated read group
			alignment.setReadGroup(settings.altReadGroupId);
		}

		//add as ready for writing
		addToContainer(_alignmentStorage, &alignment, true);

	} else if(settings.type == paired){
		//is orphan
		if(_keepOrphans){
			//add as ready for writing
			addToContainer(_alignmentStorage, &alignment, true);
		} else {
			//filter out (ignore) but write reason to bam log
			_bamFile.filterOut(alignment.name(), alignment.isReverseStrand(), alignment.readGroupId());
		}
	}
};

bool TAlignmentSplitMerger::_alignmentCanBeWrittenUnchanged(){
	return  !_recalibrate &&
			!_bamFile.curIsPaired() &&
			_alignmentStorage.empty() &&
			_rgSettings.getType(_bamFile.curReadGroupID())==unchanged;
}

//-----------------------------------------
// TOverlapQuantifier
//-----------------------------------------

void TOverlapQuantifier::quantifyOverlap(){
	//prepare counter
	coretools::TCountDistributionVector overlapDist;
	//parse BAM file
	_bamFile.startProgressReporting(1000000);
	while(_bamFile.readNextAlignment()){
		//if on new chromosome, empty storage
		if(_bamFile.chrChanged()){
			//clear storage
			_alignmentStorage.clear();
		}

		//check if first alignment in storage is too far away from current alignment
		//if yes, first alignment in storage is considered an orphan
		auto it = _alignmentStorage.begin();
		while(it != _alignmentStorage.end() && (_bamFile.curPosition() - it->alignment()) > _bamFile.maxReadLength()){
			it = _alignmentStorage.erase(it);
		}

		//check if read passed filters and is proper pair
		if(_bamFile.curPassedQC() && _bamFile.curIsProperPair()){
			//parse alignment
			BAM::TAlignment* alignment = new BAM::TAlignment;
			_bamFile.fill(*alignment);

			//check if mate is in storage.
			auto mate = findInStorage(_alignmentStorage, alignment->name());
			if(mate == _alignmentStorage.end()){
				//add alignment to storage and wait for mate
				_alignmentStorage.emplace_back(alignment, false);
			} else {
				//mate found
				if(alignment->readGroupId() != mate->alignment().readGroupId()){
					throw "Mates '" + alignment->name() + "' are in different read groups!";
				}

				//calculate overlap and fragment length and add to storage
				uint32_t overlap = _merger.determineOverlapLength(*alignment, mate->alignment()).first;
				uint32_t fragmentLength = alignment->fragmentLength();

				overlapDist.add(fragmentLength, overlap);
			}
		}

		//report
		_bamFile.printProgress();
	}

	//done parsing bam file: report
	_bamFile.printSummary(_outputName);
	_bamFile.close();

	//write distribution
	std::string filename = _outputName + "_overlapStats.txt";
	logfile().listFlush("Writing distribution of fragment length and overlap to file '" + filename + "' ...");
	const std::vector<std::string> header = {"fragmentLength", "overlap", "count"};
	coretools::TOutputFile out(filename, header);
	overlapDist.write(out);
	out.close();
	logfile().done();
};



}; //end namespace GenomeTasks

}; //end namespace AlignmentMerger
