/*
 * TAlignmentMerger.cpp
 *
 *  Created on: Mar 28, 2019
 *      Author: wegmannd
 */

#include "TBamFilter.h"

#include <math.h>
#include <stdlib.h>

#include <cstdint>
#include <exception>
#include <iostream>
#include <utility>
#include <vector>

#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TAlignment.h"
#include "TFile.h"
#include "TGenomePosition.h"
#include "TGenotypeData.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TReadGroups.h"
#include "TSequencedBase.h"
#include "counters.h"
#include "probability.h"
#include "stringFunctions.h"
#include "strongTypes.h"

namespace GenomeTasks{

using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using namespace coretools::str;

bool operator<(const uint16_t left, const TAlignmentMergerReadGroupSetting & right){
	return left < right.readGroupId;
};

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
		coretools::TInputFile in(readGroupSettingsFile, {"ReadGroup", "SeqType", "MaxCycles"});
		if(in.numCols() != 3){
			throw "Wrong number of entries in file '" + readGroupSettingsFile + "': need three columns corresponding to the read group name, read group type and max cycles!";
		}

		//read settings file
		uint16_t numNotInUse = 0;
		while(in.read(vec)){
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
						throw "Error reading file '" + in.name() + "' on line " + toString(in.lineNumber()) + ": max cycles should be a number!";
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
						throw "Error reading file '" + in.name() + "' on line " + toString(in.lineNumber()) + ": max cycles must be > 0 for read groups of type 'single'!";
					}

					//add to settings and create truncated read group
					auto it = _settings.emplace(rgId, single, maxCycles);
					it.first->altReadGroupId = readGroups.addAlternativeRG(vec[0] + "_truncated", vec[0]).id;

				} else if(vec[1] == "mixed"){
					if(maxCycles < 1){
						throw "Error reading file '" + in.name() + "' on line " + toString(in.lineNumber()) + ": max cycles must be > 0 for read groups of type 'mixed'!";
					}

					//add to settings and create truncated read group
					auto it = _settings.emplace(rgId, mixed, maxCycles);
					it.first->altReadGroupId = readGroups.addAlternativeRG(vec[0] + "_truncated", vec[0]).id;

				} else if(vec[1] == "paired"){
					_settings.emplace(rgId, paired, 0);
				} else {
					throw "Error reading file '" + in.name() + "' on line " + toString(in.lineNumber()) + ": Unknown read group type '" + vec[1] + "'! Expected 'unchanged', 'single', 'mixed' or 'paired'.";
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
// TBamFilter
//-----------------------------------------
TBamFilter::TBamFilter():TGenome_parsed(){
	//max distance between mates
	_maxDistanceBetweenMates = parameters().getParameterWithDefault<int>("acceptedDistance", 2000);
	logfile().list("Mates that are farther than " + toString(_maxDistanceBetweenMates) + " apart will be considered orphans. (parameter 'acceptedDistance')");

	//keep orphans
	if(parameters().parameterExists("keepOrphans")){
		_keepOrphans = true;
		logfile().list("Will keep keep orphaned reads. (parameter 'keepOrphans')");
	} else {
		_keepOrphans = false;
		logfile().list("Will filter out orphaned reads (use 'keepOrphans' to keep them).");
	}

	//recalibrate BAM?
	if(_genotypeLikelihoodCalculator.recalibrationChangesQualities() || parameters().parameterExists("incorporatePMD")){
		_recalibrate = true;
		logfile().list("Will write recalibrated quality scores.");
		if(parameters().parameterExists("incorporatePMD")){
			logfile().list("Probability of PMD will be reflected in new quality scores. (paramer 'incorporatePMD')");
			_incorporatePMD = true;
			if(!_genotypeLikelihoodCalculator.hasPMD()){
				throw "No PMD probabilities provided! Provide PMD probabilities or remove parameter 'incorporatePMD'.";
			}
		} else {
			_incorporatePMD = false;
			logfile().list("PMD will not be reflected in the quality scores. (recommended option. Use 'incorporatePMD' to overrule)");
		}
	} else {
		logfile().list("Will write original quality scores. (provide recalibration parameters to update quality scores)");
		_recalibrate = false;
		_incorporatePMD = false;
	}

	//set all read groups to unchanged
	_rgSettings.setAllAsUnchanged(_bamFile.readGroups());
};

void TBamFilter::_openBamFileForWriting(){
	TGenome_basic::_openBamForWriting(_outputName + "_filtered.bam", _outBam);
};

void TBamFilter::_writeAlignment(TAlignmentInStorage & it){
	//save the alignment to the bam file
	_outBam.writeAlignment(*(it->alignment));
	//delete it->alignment;
	it = _alignmentStorage.erase(it);
};

void TBamFilter::_writeOrFilterAsOrphan(TAlignmentInStorage & it){
	if(it->ready){
		_writeAlignment(it);
	} else if(_keepOrphans){
		//set as improper pair
		it->setAsNonProperPair();
		//write to BAM file
		_writeAlignment(it);
	} else {
		//write reason to bam log
		_bamFile.filterOut(it->alignment->name(), it->alignment->isSecondMate(), it->alignment->readGroupId());
		it = _alignmentStorage.erase(it);
	}
};

void TBamFilter::_writeAll(){
	//write everything and mark reads with missing mates as improper.
	//reads still in storage are no-proper pairs: write or add to black list
	TAlignmentInStorage it = _alignmentStorage.begin();
	while(it != _alignmentStorage.end()){
		_writeOrFilterAsOrphan(it);
	}
	//clear blacklist: future reads will anyways be orphans
	_blacklist.clear();
};

void TBamFilter::_writeUpTo(const genometools::TGenomePosition & position){
	//writes all that are ready or too far away
	TAlignmentInStorage it = _alignmentStorage.begin();
	while(it != _alignmentStorage.end() && (it->ready || static_cast<uint32_t>(abs(position - *it->alignment)) > _maxDistanceBetweenMates)){
		_writeOrFilterAsOrphan(it);
	}
};

BAM::TAlignment* TBamFilter::_parseIntoNewAlignment(){
	BAM::TAlignment* alignment = new BAM::TAlignment;
	_bamFile.fill(*alignment);
	if(_recalibrate){
		if(_incorporatePMD){
			alignment->parse();
			alignment->recalibrateWithPMD(_genotypeLikelihoodCalculator);
		} else {
			alignment->parse(_genotypeLikelihoodCalculator.getSequencingErrorModels());
		}
	}
	return alignment;
};

void TBamFilter::_handleMates(BAM::TAlignment* alignment, TAlignmentInStorage & mate){
	if(!alignment->isProperPair()){
		//not a proper pair: mark mate as as improper
		mate->setAsNonProperPair();
	}
	//mark both as ready for writing
	mate->ready = true;
	_alignmentStorage.emplace_back(alignment, true);
};

void TBamFilter::_handleSingle(BAM::TAlignment* alignment){
	//read is single end: add for writing
	_alignmentStorage.emplace_back(alignment, true);
};

void TBamFilter::traverseBAM(){
	//open writer
	_openBamFileForWriting();
	_bamFile.setExternalFilterReason("Orphan");

	//now parse BAM file
	_bamFile.startProgressReporting(1000000);
	while(_bamFile.readNextAlignment()){
		//if on new chromosome, empty storage
		if(_bamFile.chrChanged()){
			//write all ready currently in storage
			_writeAll();
		}

		//check if first alignment in storage is too far away from current alignment
		//if yes, first alignment in storage is considered an orphan
		_writeUpTo(_bamFile.curPosition());

		//check if read passed filters
		if(_bamFile.curPassedQC()){
			//if single end, unchanged and storage is empty: write directly
			if(!_recalibrate && !_bamFile.curIsPaired() && _alignmentStorage.empty() && _rgSettings.getType(_bamFile.curReadGroupID())==unchanged){
				_bamFile.writeCurAlignment(_outBam);
			} else {
				//parse alignment
				BAM::TAlignment* alignment = _parseIntoNewAlignment();

				//if read is paired, check for mate
				if(alignment->isPaired()){
					//if mate is in blacklist: add as improper pair for writing
					if(_blacklist.isInBlacklist(alignment->name())){
						//alignment->setIsProperPair(false);
						_alignmentStorage.emplace_back(alignment, false);
						_blacklist.remove(alignment->name());
					} else {
						//check if mate is in storage.
						auto mate = _alignmentStorage.find(alignment->name());
						if(mate == _alignmentStorage.end()){
							//no mate found
							if(alignment->isProperPair()){
								//proper pair: add to storage and wait for mate
								_alignmentStorage.emplace_back(alignment, false);
							} else {
								//improper pair: add to blacklist and ready to write
								_blacklist.add(alignment->name());
								_alignmentStorage.emplace_back(alignment, true);
							}
						} else {
							if(alignment->readGroupId() != mate->alignment->readGroupId()){
								throw "Mates '" + alignment->name() + "' are in different read groups!";
							}
							//mate found
							_handleMates(alignment, mate);
						}
					}
				} else {
					//read is single end
					_handleSingle(alignment);
				}
			}
		} else {
			//Did not pass QC: filter out
			//need to store in blacklist if is was paired
			if(_bamFile.curIsProperPair()){
				_blacklist.add(_bamFile.curName());
			}
		}

		//report
		_bamFile.printProgress();
	}

	//write reads still in storage
	_writeAll();

	//done parsing bam file: report
	_bamFile.printSummary();
	_bamFile.close();
	_outBam.close(&logfile());
};

//-----------------------------------------
// TAlignmentMergerType
//-----------------------------------------
uint16_t TAlignmentMerger::merge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	//NOTE: mate is earlier!
	//deletions and insertions are kept as is. these positions are not compared

	//check if reads overlap
	if(alignment > mate.lastAlignedPositionWithRespectToRef()){
		return 0;
	}

	//prepare (e.g. pick random number
	uint16_t numOverlap = 0;

	//go through alignments
	uint32_t fwdP = 0; uint32_t revP = 0;
	while(fwdP <= alignment.lastAlingedInternalPos() && revP <= mate.lastAlingedInternalPos()){
		
		//make sure we compare at the same position in respect to ref
		if(!alignment.isAlignedAtInternalPos(fwdP)){
			++fwdP;
		} else if(!mate.isAlignedAtInternalPos(revP)){
			++revP;
		} else if(alignment.positionInRef(fwdP) < mate.positionInRef(revP)){
			++fwdP;
		} else if(alignment.positionInRef(fwdP) > mate.positionInRef(revP)){
			++revP;
		} else {
			//at same position: merge
			++numOverlap;
			_mergeBases(alignment[fwdP], mate[revP]);

			//advance in both reads
			++fwdP; ++revP;
		}
	}

	//check if alignments changed
	if(numOverlap > 0){
		alignment.setSequenceAndQualitiesChanged();
		mate.setSequenceAndQualitiesChanged();
	}

	return numOverlap;
};

// TAlignmentMergerType_randomBase
//---------------------------------
TAlignmentMerger_randomBase::TAlignmentMerger_randomBase(const bool AdaptQuality){
	_adaptQuality = AdaptQuality;
};

void TAlignmentMerger_randomBase::_mergeBasesCore(BAM::TSequencedBase & bestBase, BAM::TSequencedBase & worstBase){
	if(_adaptQuality){
		auto likelihood = GenotypeLikelihoods::fromError(
		    bestBase.base, (coretools::Probability)bestBase.recalibratedQualityAsPhredInt);
		const auto worst = GenotypeLikelihoods::fromError(
		    worstBase.base, (coretools::Probability)worstBase.recalibratedQualityAsPhredInt);

		for (auto b = genometools::Base::min; b < genometools::Base::max; ++b) likelihood[b] *= worst[b];

		coretools::normalize(likelihood);
		bestBase.recalibratedQualityAsPhredInt = likelihood[bestBase.base].complement();
	}

	//set other to missing
	worstBase.recalibratedQualityAsPhredInt = 0.0;
	worstBase.base = genometools::Base::N;
};

void TAlignmentMerger_randomBase::_mergeBases(BAM::TSequencedBase & alignment, BAM::TSequencedBase & mate){
	if(randomGenerator().pickOneOfTwo()){
		_mergeBasesCore(mate, alignment);
	} else {
		_mergeBasesCore(alignment, mate);
	}
};

// TAlignmentMergerType_randomRead
//---------------------------------
TAlignmentMerger_randomRead::TAlignmentMerger_randomRead(const bool AdaptQuality):TAlignmentMerger_randomBase(AdaptQuality){
	_keepMate = false;
};

void TAlignmentMerger_randomRead::_mergeBases(BAM::TSequencedBase & alignment, BAM::TSequencedBase & mate){
	if(_keepMate){
		_mergeBasesCore(mate, alignment);
	} else {
		_mergeBasesCore(alignment, mate);
	}
};

uint16_t TAlignmentMerger_randomRead::merge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	_keepMate = randomGenerator().pickOneOfTwo();
	return TAlignmentMerger::merge(alignment, mate);
};

// TAlignmentMergerType_highestQuality
//---------------------------------
TAlignmentMerger_highestQuality::TAlignmentMerger_highestQuality(const bool AdaptQuality):TAlignmentMerger_randomBase(AdaptQuality){};

void TAlignmentMerger_highestQuality::_mergeBases(BAM::TSequencedBase & alignment, BAM::TSequencedBase & mate){
	if(mate.recalibratedQualityAsPhredInt > alignment.recalibratedQualityAsPhredInt){
		_mergeBasesCore(mate, alignment);
	} else if(alignment.recalibratedQualityAsPhredInt > mate.recalibratedQualityAsPhredInt){
		_mergeBasesCore(alignment, mate);
	} else {
		//pick randomly
		TAlignmentMerger_randomBase::_mergeBases(alignment, mate);
	}
};

//-----------------------------------------
// TAlignmentSplitMerger
//-----------------------------------------
TAlignmentSplitMerger::TAlignmentSplitMerger() : TBamFilter(){
	//parse read group settings
	_rgSettings.initialize(_bamFile.readGroupsMutable());

	//other settings
	_numReadsMerged = 0;
	_numBasesMerged = 0;

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
	*/

	//decide if we update quality score
	bool adaptQuality;
	if(parameters().parameterExists("updateQuality")){
		adaptQuality = true;
		logfile().list("Will update quality scores of preferred bases to reflect information from overlapping bases.");
	} else {
		adaptQuality = false;
		logfile().list("Will keep original quality scores of the preferred bases (use updateQuality to update quality scores).");
	}

	//set merging method
	//TODO: update wiki to reflect change in names
	std::string method = parameters().getParameterWithDefault<std::string>("mergingMethod", "randomRead");
	if(method == "none"){
		_merger = std::make_unique<TAlignmentMerger>();
		logfile().list("Merging method: no merging.");
	} else if (method == "randomRead"){
		_merger = std::make_unique<TAlignmentMerger_randomRead>(adaptQuality);
		logfile().list("Merging method: will keep random read for all overlapping positions");
	} else if(method == "randomBase"){
		_merger = std::make_unique<TAlignmentMerger_randomBase>(adaptQuality);
		logfile().list("Merging method: will keep random base at each overlapping position.");
	} else if(method == "highestQuality"){
		_merger = std::make_unique<TAlignmentMerger_highestQuality>(adaptQuality);
		logfile().list("Merging method: will keep base with highest quality at overlapping positions.");
	} else {
		throw "Unknown merging method " + method + "! Use 'none', 'randomRead', 'randomBase' and 'highestQuality'.";
	}
};

void TAlignmentSplitMerger::_openBamFileForWriting(){
	TGenome_basic::_openBamForWriting(_outputName + "_splitMerged.bam", _outBam);
};

void TAlignmentSplitMerger::_handleMates(BAM::TAlignment* alignment, TAlignmentInStorage & mate){
	ReadGroupType type = _rgSettings.getType(alignment->readGroupId());

	if(type == single){
		throw "Paired reads found in single-end read group '" + _bamFile.readGroups().getName(alignment->readGroupId()) + "'! Is this a 'mixed' read group?";
	} else if(!alignment->isProperPair()){
		//not a proper pair: mark mate as as improper too
		mate->setAsNonProperPair();
	} else if(type == paired || type == mixed){
		//attempt merging: make sure alignments are parsed
		//Note: if we recalibrate, they were already parsed
		if(!alignment->isParsed()){
			alignment->parse();
		}
		if(!mate->alignment->isParsed()){
			mate->alignment->parse();
		}

		uint16_t overlap = _merger->merge(*alignment, *mate->alignment);
		if(overlap > 0){
			++_numReadsMerged;
			_numBasesMerged += overlap;
		}
	}

	//mark both as ready for writing
	mate->ready = true;
	_alignmentStorage.emplace_back(alignment, true);
};

void TAlignmentSplitMerger::_handleSingle(BAM::TAlignment* alignment){
	const TAlignmentMergerReadGroupSetting& settings = _rgSettings.getSettings(alignment->readGroupId());

	if(settings.type == unchanged){
		//add as ready for writing
		_alignmentStorage.emplace_back(alignment, true);
	} else if(settings.type == single || settings.type == mixed){
		//truncate
		if(alignment->parsedLength() > settings.maxCycles && !_allowForLarger){
			throw("Length of read " + alignment->name() + " is > max cycles for its read group (" + toString(settings.maxCycles) + ")! Use parameter 'allowForLarger' to ignore and put read in truncated read group.");
		} else if(alignment->parsedLength() >= settings.maxCycles - _considerAtMaxUpToDist){
			//add to truncated read group
			alignment->setReadGroup(settings.altReadGroupId);
		}

		//add as ready for writing
		_alignmentStorage.emplace_back(alignment, true);

	} else if(settings.type == paired){
		//is orphan
		if(_keepOrphans){
			//add as ready for writing
			_alignmentStorage.emplace_back(alignment, true);
		} else {
			//filter out (ignore) but write reason to bam log
			_bamFile.filterOut(alignment->name(), alignment->isReverseStrand(), alignment->readGroupId());
		}
	}
};

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
		TAlignmentInStorage it = _alignmentStorage.begin();
		while(it != _alignmentStorage.end() && (_bamFile.curPosition() - *it->alignment) > _bamFile.maxReadLength()){
			it = _alignmentStorage.erase(it);
		}

		//check if read passed filters and is proper pair
		if(_bamFile.curPassedQC() && _bamFile.curIsProperPair()){
			//parse alignment
			BAM::TAlignment* alignment = new BAM::TAlignment;
			_bamFile.fill(*alignment);

			//check if mate is in storage.
			auto mate = _alignmentStorage.find(alignment->name());
			if(mate == _alignmentStorage.end()){
				//add alignment to storage and wait for mate
				_alignmentStorage.emplace_back(alignment, false);
			} else {
				//mate found
				if(alignment->readGroupId() != mate->alignment->readGroupId()){
					throw "Mates '" + alignment->name() + "' are in different read groups!";
				}

				//calculate overlap and fragment length and add to storage
				uint16_t overlap = _merger.merge(*alignment, *mate->alignment);
				uint16_t fragmentLength = alignment->parsedLength() + mate->alignment->parsedLength() - overlap;

				overlapDist.add(fragmentLength, overlap);
			}
		}

		//report
		_bamFile.printProgress();
	}

	//done parsing bam file: report
	_bamFile.printSummary();
	_bamFile.close();

	//write distribution
	std::string filename = _outputName + "_overlapStats.txt";
	logfile().listFlush("Writing distribution of fragment length and overlap to file '" + filename + "' ...");
	const std::vector<std::string> header = {"fragmentLength", "overlap"};
	coretools::TOutputFile out(filename, header);
	overlapDist.write(out);
	out.close();
	logfile().done();
};



}; //end namespace





