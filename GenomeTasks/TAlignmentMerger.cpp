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

#include "coretools/Containers/TStrongArray.h"
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


namespace impl{
	void callMergeFunction(BAM::TAlignment & alignment, BAM::TAlignment & mate, size_t overlapLength){
		size_t lengthMappedBeforeMerge = alignment.cigar().lengthMapped();
		size_t mappedBasesClipped = 0;
		//the reads are merged by softclipping bases, for which a new cigar string is constructed
		alignment.merge(overlapLength, mappedBasesClipped);
		if(alignment.isReverseStrand()){
			//if the reverse strand gets merged, its position in the BAM-file as well as the position of the mate of the first read need to be adjusted to account for the length of the added softclips on left side
			if (overlapLength < lengthMappedBeforeMerge) {
				alignment.moveOnRef(alignment.position() + mappedBasesClipped);
				mate.setMateGenomicPosition(alignment);
			} else {
				//if the overlap is larger than the read being merged, the position can only be moved up by the amount of aligned bases in the read
				alignment.moveOnRef(alignment.position() + lengthMappedBeforeMerge);
				mate.setMateGenomicPosition(alignment);
			}
		}
	}

	std::pair<genometools::PhredIntProbability,genometools::PhredIntProbability> minQual(const BAM::TAlignment & firstRead, const BAM::TAlignment & secondRead) {
		//base iterator starts at first position of the reverse strand, then increments until it reaches either the last aligned position of itself or the forward read
		std::vector<BAM::TSequencedBase>::const_iterator baseIterator = secondRead.begin();
		size_t internalPos = 0;
		//use the recalibrated quality
		genometools::PhredIntProbability secondReadMinQual = baseIterator->recalibratedQualityAsPhredInt;

		while(!secondRead.isAlignedAtInternalPos(internalPos) || (secondRead.positionInRef(internalPos).position() != firstRead.lastAlignedPositionWithRespectToRef().position() && secondRead.positionInRef(internalPos).position() != secondRead.lastAlignedPositionWithRespectToRef().position())){
			//if the base at the current position of the iterator is aligned and its PhredIntProbability is higher (therefore the error-probability is higher and the quality is lower)
			//save this number as the new minimum quality
			if(secondRead.isAlignedAtInternalPos(internalPos)){
				if (baseIterator->recalibratedQualityAsPhredInt > secondReadMinQual)
					secondReadMinQual = baseIterator->recalibratedQualityAsPhredInt;
			}
			baseIterator++;
			internalPos++;
		}
		//base iterator starts at last position of the forward strand, then decrements until it reaches either the first aligned position of itself or the forward read
		std::vector<BAM::TSequencedBase>::const_reverse_iterator baseIteratorReverse = firstRead.rbegin();
		internalPos = firstRead.getLastInternalPos();
		genometools::PhredIntProbability firstReadMinQual = baseIteratorReverse->recalibratedQualityAsPhredInt;
		while (!firstRead.isAlignedAtInternalPos(internalPos) || (firstRead.positionInRef(internalPos).position() != secondRead.position() && firstRead.positionInRef(internalPos).position() != firstRead.position())) {
			if(firstRead.isAlignedAtInternalPos(internalPos)){
				if (baseIteratorReverse->recalibratedQualityAsPhredInt > firstReadMinQual)
					firstReadMinQual = baseIteratorReverse->recalibratedQualityAsPhredInt;
			}
			baseIteratorReverse++;
			internalPos--;
		}
		//returns a pair of the minimum qualities of both reads (in the overlap)
		return std::make_pair(firstReadMinQual,secondReadMinQual);
	}

	std::pair<genometools::PhredIntProbability,genometools::PhredIntProbability> getMinQuals(const BAM::TAlignment & alignment, const BAM::TAlignment & mate) {
		//this if-statement ensures that the order of the pair of minimum qualities that is returned matches the order in which the two reads were passed to the function
		//this is necessary because the minQual function needs the two reads in the order of (forwardStrand, reverseStrand)
		if (!alignment.isReverseStrand()){
			return minQual(alignment,mate);
		} else {
			std::pair<genometools::PhredIntProbability,genometools::PhredIntProbability> flippedResult = minQual(mate,alignment);
			return std::make_pair(flippedResult.second,flippedResult.first);
		}
	}

	genometools::PhredIntProbability determineQualAtSingleBase(const BAM::TAlignment & alignment, size_t numBasesFromEnd, bool isFirst){
		//for the first read, we start at the last base of the read and go to the base at the center of the overlap and return its quality-value
		if (isFirst){
			std::vector<BAM::TSequencedBase>::const_reverse_iterator baseIterator = alignment.rbegin() + numBasesFromEnd;
			return baseIterator->recalibratedQualityAsPhredInt;
		} else {
			//for the second read, we start at the first base of the read
			std::vector<BAM::TSequencedBase>::const_iterator baseIterator = alignment.begin() + numBasesFromEnd;
			return baseIterator->recalibratedQualityAsPhredInt;
		}
	}

	void compareQualities(const genometools::PhredIntProbability & alignmentQual, const genometools::PhredIntProbability & mateQual, size_t & alignmentOverlapLength, size_t & mateOverlapLength){
		if (alignmentQual < mateQual) 
			mateOverlapLength++;
		else if (alignmentQual > mateQual) 
			alignmentOverlapLength++;
		else { 
			//if both qualities are the same, pick one randomly
			if (randomGenerator().pickOneOfTwo()) 
				mateOverlapLength++;
			else 
				alignmentOverlapLength++;
		}
	}

	void checkIfReverseStrandFirst(const BAM::TAlignment & alignment, const BAM::TAlignment & mate, size_t & alignmentOverlapLength, size_t & mateOverlapLength, bool alignmentIsFirst){
		//if the reverse strand is the first read
		//we need to add the distance between the start positions to the overlap length of the reverse strand
		//and the distance between the end positions to the overlap length of the forward strand
		if (alignmentIsFirst && alignment.isReverseStrand()){
			alignmentOverlapLength += mate.position() - alignment.position();
			mateOverlapLength += (mate.position() + mate.cigar().lengthMapped()) - (alignment.position() + alignment.cigar().lengthMapped());
		} else if (!alignmentIsFirst && mate.isReverseStrand()) {
			mateOverlapLength += alignment.position() - mate.position();
			alignmentOverlapLength += (alignment.position() + alignment.cigar().lengthMapped()) - (mate.position() + mate.cigar().lengthMapped());
		}
	}

	void mergeOddOverlap(BAM::TAlignment & firstRead, BAM::TAlignment & secondRead, size_t halfOverlap){
		//calculate qualities for the position at the center of the overlap for both reads
		genometools::PhredIntProbability firstReadQual = determineQualAtSingleBase(firstRead, halfOverlap, true);
		genometools::PhredIntProbability secondReadQual = determineQualAtSingleBase(secondRead, halfOverlap, false);

		size_t firstReadOverlapLength = halfOverlap;
		size_t secondReadOverlapLength = halfOverlap;
		//compare qualities and add 1 to the overlap length of the read with the lower quality base at the center of the overlap
		compareQualities(firstReadQual, secondReadQual, firstReadOverlapLength, secondReadOverlapLength);

		//check if the reverse strand comes first and adjust the overlap lengths accordingly if that is the case
		checkIfReverseStrandFirst(firstRead, secondRead, firstReadOverlapLength, secondReadOverlapLength, true);
		//merge
		callMergeFunction(firstRead, secondRead, firstReadOverlapLength);
		callMergeFunction(secondRead, firstRead, secondReadOverlapLength);
	}

	void mergeOverlapLargerThanRead(BAM::TAlignment & largerRead, BAM::TAlignment & smallerRead, size_t largerReadOverlapLength, size_t smallerReadOverlapLength){
		//if the overlap is divisible by two, you can just merge
		if (smallerRead.cigar().lengthMapped() % 2 == 0){
			callMergeFunction(smallerRead, largerRead, smallerReadOverlapLength);
			callMergeFunction(largerRead, smallerRead, largerReadOverlapLength);
		} else {
			//if the overlap is not divisible by two, we first need to determine the quality at the position in the center of the overlap for both reads
			genometools::PhredIntProbability smallerReadQual = determineQualAtSingleBase(smallerRead, smallerReadOverlapLength, !smallerRead.isReverseStrand());
			genometools::PhredIntProbability largerReadQual = determineQualAtSingleBase(largerRead, largerReadOverlapLength, !largerRead.isReverseStrand());
		
			compareQualities(smallerReadQual, largerReadQual, smallerReadOverlapLength, largerReadOverlapLength);
			callMergeFunction(smallerRead, largerRead, smallerReadOverlapLength);
			callMergeFunction(largerRead, smallerRead, largerReadOverlapLength);
		}
	}

	void handleOverlapLargerThanRead(BAM::TAlignment & largerRead, BAM::TAlignment & smallerRead){
		//since the larger read fully eclipses the smaller read, the overlap is now as large as the number of aligned positions in the smaller read
		size_t halfOverlap = smallerRead.cigar().lengthMapped() / 2;
		size_t smallerReadOverlapLength = halfOverlap;
		size_t largerReadOverlapLength = halfOverlap;
		//if the larger read is the forward strand, add the distance between the ends of both reads to its overlap
		//if it is the reverse read, add the distance between both starts
		if (!largerRead.isReverseStrand()){
			largerReadOverlapLength += (largerRead.position() + largerRead.cigar().lengthMapped()) - (smallerRead.position() + smallerRead.cigar().lengthMapped());
		} else {
			largerReadOverlapLength += smallerRead.position() - largerRead.position();
		}

		//if both reads are forward/reverse, the smaller read needs to be set to reverse/forward so different sides of the reads are soft-clipped
		if ((largerRead.isReverseStrand() && smallerRead.isReverseStrand())){
			smallerRead.setIsReverseStrand(false);
			mergeOverlapLargerThanRead(largerRead, smallerRead, largerReadOverlapLength, smallerReadOverlapLength);
			smallerRead.setIsReverseStrand(true);
		} else if ((!largerRead.isReverseStrand() && !smallerRead.isReverseStrand())) {
			smallerRead.setIsReverseStrand(true);
			mergeOverlapLargerThanRead(largerRead, smallerRead, largerReadOverlapLength, smallerReadOverlapLength);
			smallerRead.setIsReverseStrand(false);
		} else {
			mergeOverlapLargerThanRead(largerRead, smallerRead, largerReadOverlapLength, smallerReadOverlapLength);
		}
	}
}

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
			for (size_t rg = 0; rg < readGroups.size(); ++rg) { _settings.emplace(rg, ReadGroupType::paired, 0); }
		} else {
			    // will merge a subset and treat others as unchanged
			    std::vector<std::string> vec;
			    fillContainerFromString(pairedRG, vec, ',');
			    logfile().listFlush("Parsing read group names from parameter 'pairedReadGroups' ...");

			    // get IDs
			    std::set<size_t, std::less<>> pairedIds;
			    for (auto n : vec) { pairedIds.emplace(readGroups.getId(n)); }

			    // add as unchanged or paired
			    for (size_t rg = 0; rg < readGroups.size(); ++rg) {
				    if (pairedIds.find(rg) == pairedIds.end()) {
					    _settings.emplace(rg, ReadGroupType::unchanged, 0);
				    } else {
					    _settings.emplace(rg, ReadGroupType::paired, 0);
				    }
			    }
			    logfile().done();
		}
		_printSummary();
	} else {
		//do we have to ignore read groups present in file?
		std::vector<std::string> vec;
		std::set<size_t> readGroupsToIgnore;
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
		size_t numNotInUse = 0;
		while(in.read(vec)){
			//ignore "allReadGroups" from BAMD output
			if (vec[0] != "allReadGroups"){
				//get read group ID
				//read groups not in use will get a warning and be ignored
				size_t rgId = readGroups.getId(vec[0]);
				if(!readGroups.readGroupInUse(rgId)){
					++numNotInUse;
				} else {
					//parse max cycles
					size_t maxCycles = 0;
					if(vec[2] != "NA" && vec[2] != "-"){
						if(!stringContainsOnlyNumbers(vec[2])){
							UERROR("Error reading file '", in.name(), "' on line ", in.curLine(), ": max cycles should be a number!");
						}
						maxCycles = fromString<int>(vec[2]);
					}

					//check for duplicate entries
					if(_settings.find(rgId) != _settings.end()){
						UERROR("Duplicate entry in file '", in.name(), "' for read group '", vec[0], "'!");
					}

					//act based on seqeuncing type (second column). Ignored read groups will be marked as "unchanged"
					if(readGroupsToIgnore.find(rgId) != readGroupsToIgnore.end() || vec[1] == "unchanged"){
						_settings.emplace(rgId, ReadGroupType::unchanged, 0);
					} else if(vec[1] == "single"){
						if(maxCycles < 1){
							UERROR("Error reading file '", in.name(), "' on line ", in.curLine(), ": max cycles must be > 0 for read groups of type 'single'!");
						}

						//add to settings and create truncated read group
						_settings.emplace(rgId, readGroups.addAlternativeRG(vec[0] + "_truncated", vec[0]).id(), ReadGroupType::single, maxCycles);

					} else if(vec[1] == "mixed"){
						if(maxCycles < 1){
							UERROR("Error reading file '", in.name(), "' on line ", in.curLine(), ": max cycles must be > 0 for read groups of type 'mixed'!");
						}

						//add to settings and create truncated read group
						_settings.emplace(rgId, readGroups.addAlternativeRG(vec[0] + "_truncated", vec[0]).id(), ReadGroupType::mixed, maxCycles);

					} else if(vec[1] == "paired"){
						_settings.emplace(rgId, ReadGroupType::paired, 0);
					} else {
						UERROR("Error reading file '", in.name(), "' on line ", in.curLine(), ": Unknown read group type '", vec[1], "'! Expected 'unchanged', 'single', 'mixed' or 'paired'.");
					}
				}
			}
		}

		//set missing read groups to "unchanged"
		for(size_t rg=0; rg < readGroups.size(); ++rg){
			if(readGroups.readGroupInUse(rg) && _settings.find(rg)==_settings.end()){
				_settings.emplace(rg, ReadGroupType::unchanged, 0);
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
	coretools::TStrongArray<size_t, ReadGroupType> counts;
	for(auto& s : _settings){
		++counts[s.type];
	}

	//summarize
	if(counts[ReadGroupType::unchanged] > 0){ logfile().conclude(counts[ReadGroupType::unchanged], " read groups will remain unchanged."); }
	if(counts[ReadGroupType::single] > 0   ){ logfile().conclude(counts[ReadGroupType::single], " single-end read groups will be split."); }
	if(counts[ReadGroupType::mixed] > 0    ){ logfile().conclude(counts[ReadGroupType::mixed], " mixed read groups will be split and merged."); }
	if(counts[ReadGroupType::paired] > 0   ){ logfile().conclude(counts[ReadGroupType::paired], " paired read groups to be merged."); }
};

void TAlignmentMergerReadGroupSettings::setAllAsUnchanged(const BAM::TReadGroups & readGroups){
	_settings.clear();
	for(size_t rg=0; rg < readGroups.size(); ++rg){
		if(readGroups.readGroupInUse(rg)){
			_settings.emplace(rg, ReadGroupType::unchanged, 0);
		}
	}
};

bool TAlignmentMergerReadGroupSettings::needTruncation() const{
	for(auto& s : _settings){
		if(s.type == ReadGroupType::single || s.type == ReadGroupType::mixed)
			return true;
	}
	return false;
};


bool TAlignmentMergerReadGroupSettings::needsMerging() const{
	for(const auto& s : _settings){
		if(s.type == ReadGroupType::paired || s.type == ReadGroupType::mixed)
			return true;
	}
	return false;
};

ReadGroupType TAlignmentMergerReadGroupSettings::getType(size_t readGroupId) const{
	return _settings.find(readGroupId)->type;
};

size_t TAlignmentMergerReadGroupSettings::getMaxCycles(size_t readGroupId) const{
	return _settings.find(readGroupId)->maxCycles;
};

const TAlignmentMergerReadGroupSetting& TAlignmentMergerReadGroupSettings::getSettings(size_t readGroupId) const{
	return *_settings.find(readGroupId);
};

//-----------------------------------------
// TAlignmentMergerType
//-----------------------------------------
size_t TAlignmentMerger::merge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	bool alignmentFirst = mate > alignment;
	size_t overlapLength = 0;
	if ((alignment.isReverseStrand() && mate.isReverseStrand())){
		alignmentFirst ? alignment.setIsReverseStrand(false) : mate.setIsReverseStrand(false);
		overlapLength = overlapLengthAndMerge(alignment, mate);
		!alignment.isReverseStrand() ? alignment.setIsReverseStrand(true) : mate.setIsReverseStrand(true);
	} else if ((!alignment.isReverseStrand() && !mate.isReverseStrand())) {
		alignmentFirst ? mate.setIsReverseStrand(true) : alignment.setIsReverseStrand(true);
		overlapLength = overlapLengthAndMerge(alignment, mate);
		alignment.isReverseStrand() ? alignment.setIsReverseStrand(false) : mate.setIsReverseStrand(false);
	} else {
		overlapLength = overlapLengthAndMerge(alignment, mate);
	}
	return overlapLength;
};

size_t TAlignmentMerger::determineOverlapLength(const BAM::TAlignment & alignment, const BAM::TAlignment & mate){
	//No overlap if either one of the reads is already completely softclipped
	if (alignment.cigar().lengthAligned() == 0 || mate.cigar().lengthAligned() == 0)
		return 0;
	//if-statements check for overlap using the first aligned position in the BAM-file and the aligned length given in the cigar string
	if (alignment.isReverseStrand()){
		if (alignment.position() < mate.position() + mate.cigar().lengthMapped()){
			size_t overlapLength = mate.position() + mate.cigar().lengthMapped() - alignment.position();
			return overlapLength;
		} return 0;
	} else {
		if (mate.position() < alignment.position() + alignment.cigar().lengthMapped()){
			size_t overlapLength = alignment.position() + alignment.cigar().lengthMapped() - mate.position();
			return overlapLength;
		} return 0;
	}
}

size_t TAlignmentMerger::overlapLengthAndMerge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	//check if reads overlap
	size_t overlapLength = determineOverlapLength(alignment, mate);
	//if they do -> merge
	if (overlapLength > 0){
		impl::callMergeFunction(alignment, mate, overlapLength);
	}
	return overlapLength;
}

// TAlignmentMergerType_randomRead
//---------------------------------
TAlignmentMerger_randomRead::TAlignmentMerger_randomRead():TAlignmentMerger(){};

size_t TAlignmentMerger_randomRead::merge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	//randomly either merge mate or alignment
	if (randomGenerator().pickOneOfTwo())
		return TAlignmentMerger::merge(mate, alignment);
	//else
	return TAlignmentMerger::merge(alignment, mate);
};

// TAlignmentMergerType_middle
//---------------------------------
TAlignmentMerger_middle::TAlignmentMerger_middle():TAlignmentMerger(){};

size_t TAlignmentMerger_middle::merge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	const std::pair<size_t,bool> overlapLength = determineOverlapLength(alignment, mate);
	if (overlapLength.first > 0){
		//in case one of the reads fully covers its mate 
		if (overlapLength.first >= alignment.cigar().lengthMapped()) {
			impl::handleOverlapLargerThanRead(mate, alignment);
			return overlapLength.first;
		} else if (overlapLength.first >= mate.cigar().lengthMapped()){
			impl::handleOverlapLargerThanRead(alignment, mate);
			return overlapLength.first;
		}
		//if both reads are forward/reverse strands, we need to change one of their directions before merging so we don't clip the same side twice
		if ((alignment.isReverseStrand() && mate.isReverseStrand())){
			!overlapLength.second ? mate.setIsReverseStrand(false) : alignment.setIsReverseStrand(false);
			sameDirectionMerge(alignment, mate, overlapLength);
			!alignment.isReverseStrand() ? alignment.setIsReverseStrand(true) : mate.setIsReverseStrand(true);
		} else if ((!alignment.isReverseStrand() && !mate.isReverseStrand())) {
			!overlapLength.second ? alignment.setIsReverseStrand(true) : mate.setIsReverseStrand(true);
			sameDirectionMerge(alignment, mate, overlapLength);
			alignment.isReverseStrand() ? alignment.setIsReverseStrand(false) : mate.setIsReverseStrand(false);
		} else {
			size_t halfOverlap = overlapLength.first / 2;
			//if the overlap length is divisible by 2, we only need to check if the reverse strand is first before merging
			if (overlapLength.first % 2 == 0){
				size_t alignmentOverlapLength = halfOverlap;
				size_t mateOverlapLength = halfOverlap;
				impl::checkIfReverseStrandFirst(alignment, mate, alignmentOverlapLength, mateOverlapLength, overlapLength.second);

				impl::callMergeFunction(alignment, mate, alignmentOverlapLength);
				impl::callMergeFunction(mate, alignment, mateOverlapLength);
			} else {
				//this if-statement is necessary because the mergeOddOverlap-function needs the alignments in the order of (firstRead, secondRead) 
				//to calculate the qualities at the center position
				if (overlapLength.second == true) 
					impl::mergeOddOverlap(alignment, mate, halfOverlap);
				else
					impl::mergeOddOverlap(mate, alignment, halfOverlap);
			}
		}
	} 
	return overlapLength.first;
};

std::pair<size_t,bool> TAlignmentMerger_middle::determineOverlapLength(const BAM::TAlignment & alignment, const BAM::TAlignment & mate){
	//this function works like the one it overrides, except that it checks which read comes first instead of which is the reverse strand
	//also it returns a boolean which states whether the alignment is the first read
	if (alignment.cigar().lengthAligned() == 0 || mate.cigar().lengthAligned() == 0)
		return std::make_pair(0,true);
	if (alignment > mate){
		if (alignment.position() < mate.position() + mate.cigar().lengthMapped()){
			size_t overlapLength = mate.position() + mate.cigar().lengthMapped() - alignment.position();
			bool alignmentFirst = false;
			return std::make_pair(overlapLength,alignmentFirst);
		} 
		//else
			return std::make_pair(0,false);
	} else if (mate > alignment) {
		if (mate.position() < alignment.position() + alignment.cigar().lengthMapped()){
			size_t overlapLength = alignment.position() + alignment.cigar().lengthMapped() - mate.position();
			bool alignmentFirst = true;
			return std::make_pair(overlapLength,alignmentFirst);
		} 
		//else
			return std::make_pair(0,true);
	} else {
		//if both reads have the same start position, the overlap is the lower mapped length of the two
		size_t overlapLength;
		(alignment.cigar().lengthMapped() > mate.cigar().lengthMapped()) ? overlapLength = mate.cigar().lengthMapped() : overlapLength = alignment.cigar().lengthMapped();
		bool alignmentFirst = randomGenerator().pickOneOfTwo();
		return std::make_pair(overlapLength,alignmentFirst);
	}
}

void TAlignmentMerger_middle::sameDirectionMerge(BAM::TAlignment & alignment, BAM::TAlignment & mate, std::pair<size_t,bool> overlapLength){
size_t halfOverlap = overlapLength.first / 2;
if (overlapLength.first % 2 == 0){
	size_t alignmentOverlapLength = halfOverlap;
	size_t mateOverlapLength = halfOverlap;

	impl::callMergeFunction(alignment, mate, alignmentOverlapLength);
	impl::callMergeFunction(mate, alignment, mateOverlapLength);
} else {
	//this if-statement is necessary because the mergeOddOverlap-function needs the alignments in the order of (firstRead, secondRead) 
	//to calculate the qualities at the center position
	if (overlapLength.second == true) 
		impl::mergeOddOverlap(alignment, mate, halfOverlap);
	else
		impl::mergeOddOverlap(mate, alignment, halfOverlap);
	}
}


// TAlignmentMergerType_firstMate
//---------------------------------
TAlignmentMerger_firstMate::TAlignmentMerger_firstMate():TAlignmentMerger(){};

size_t TAlignmentMerger_firstMate::merge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	//always merge first mate, keep second mate
	if (alignment.isSecondMate())
		return TAlignmentMerger::merge(mate, alignment);
	else
		return TAlignmentMerger::merge(alignment, mate);
};

// TAlignmentMergerType_secondMate
//---------------------------------
TAlignmentMerger_secondMate::TAlignmentMerger_secondMate():TAlignmentMerger(){};

size_t TAlignmentMerger_secondMate::merge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	//always merge second mate, keep first mate
	if (alignment.isSecondMate())
		return TAlignmentMerger::merge(alignment, mate);
	else
		return TAlignmentMerger::merge(mate, alignment);
};

// TAlignmentMergerType_highestQuality
//---------------------------------
TAlignmentMerger_highestQuality::TAlignmentMerger_highestQuality():TAlignmentMerger(){};

size_t TAlignmentMerger_highestQuality::merge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
	bool alignmentFirst = mate > alignment;
	size_t overlapLength = 0;
	if ((alignment.isReverseStrand() && mate.isReverseStrand())){
		alignmentFirst ? alignment.setIsReverseStrand(false) : mate.setIsReverseStrand(false);
		overlapLength = overlapLengthAndMerge(alignment, mate);
		!alignment.isReverseStrand() ? alignment.setIsReverseStrand(true) : mate.setIsReverseStrand(true);
	} else if ((!alignment.isReverseStrand() && !mate.isReverseStrand())) {
		alignmentFirst ? mate.setIsReverseStrand(true) : alignment.setIsReverseStrand(true);
		overlapLength = overlapLengthAndMerge(alignment, mate);
		alignment.isReverseStrand() ? alignment.setIsReverseStrand(false) : mate.setIsReverseStrand(false);
	} else {
		overlapLength = overlapLengthAndMerge(alignment, mate);
	}
	return overlapLength;
}

size_t TAlignmentMerger_highestQuality::overlapLengthAndMerge(BAM::TAlignment & alignment, BAM::TAlignment & mate){
		//check if reads overlap
	size_t overlapLength = TAlignmentMerger::determineOverlapLength(alignment, mate);
	//if they do -> calculate minimum quality of each read in the overlap
	if (overlapLength > 0) {
		std::pair<genometools::PhredIntProbability,genometools::PhredIntProbability> minQuals = impl::getMinQuals(alignment, mate);
		//merge the read with the lower minimum quality (the read with the higher PhredIntProbability has the higher error-rate and the lower quality)
		if (minQuals.first < minQuals.second){
			impl::callMergeFunction(mate, alignment, overlapLength);
		} else if (minQuals.second < minQuals.first){
			impl::callMergeFunction(alignment, mate, overlapLength);
		} else {
				//if both reads have equal minimum qualities, randomly choose one
				if (randomGenerator().pickOneOfTwo()){
					impl::callMergeFunction(mate, alignment, overlapLength);
				} else {
					impl::callMergeFunction(alignment, mate, overlapLength);
				}
		}
	} 
	return overlapLength;
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

	//set merging method
	//TODO: update wiki to reflect change in names
	std::string method = parameters().getParameterWithDefault<std::string>("mergingMethod", "middle");
	if(method == "none"){
		_merger = std::make_unique<TAlignmentMerger>();
		logfile().list("Merging method: no merging. (parameter 'mergingMethod')");
	} else if (method == "randomRead"){
		_merger = std::make_unique<TAlignmentMerger_randomRead>();
		logfile().list("Merging method: will keep random read for all overlapping positions. (parameter 'mergingMethod')");
	} else if(method == "highestQuality"){
		_merger = std::make_unique<TAlignmentMerger_highestQuality>();
		logfile().list("Merging method: will keep read with highest minimum quality at overlapping positions. (parameter 'mergingMethod')");
	} else if(method == "firstMate"){
		_merger = std::make_unique<TAlignmentMerger_firstMate>();
		logfile().list("Merging method: will merge first mate at overlapping positions. (parameter 'mergingMethod')");
	} else if(method == "secondMate"){
		_merger = std::make_unique<TAlignmentMerger_secondMate>();
		logfile().list("Merging method: will merge second mate at overlapping positions. (parameter 'mergingMethod')");
	} else if (method == "middle"){
		_merger = std::make_unique<TAlignmentMerger_middle>();
		logfile().list("Merging method: will keep half of the overlapping positions of each mate. (parameter 'mergingMethod')");
	} else {
		UERROR("Unknown merging method ", method, "! Use 'none', 'middle', 'firstMate', 'secondMate', 'randomRead' or 'highestQuality'.");
	}
};

void TAlignmentSplitMerger::_openBamFileForWriting(){
	TGenome_basic::_openBamForWriting(_outputName + "_splitMerged.bam", _outBam);
};

void TAlignmentSplitMerger::_handleMates(BAM::TAlignment & alignment, TAlignmentStorageSortedIterator mate){
	ReadGroupType type = _rgSettings.getType(alignment.readGroupId());

	if(type == ReadGroupType::single){
		UERROR("Paired reads found in single-end read group '", _bamFile.readGroups().getName(alignment.readGroupId()), "'! Is this a 'mixed' read group?");
	} else if(!alignment.isProperPair()){
		//not a proper pair: mark mate as as improper too
		mate->setAsNonProperPair();
		mate->makeReady();
	} else if(type == ReadGroupType::paired || type == ReadGroupType::mixed){
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

	if(settings.type == ReadGroupType::unchanged){
		//add as ready for writing
		addToContainer(_alignmentStorage, &alignment, true);
	} else if(settings.type == ReadGroupType::single || settings.type == ReadGroupType::mixed){
		//truncate
		if(!_allowForLarger && alignment.length() > settings.maxCycles){
			UERROR("Length of read ", alignment.name(), " is > max cycles for its read group (",settings.maxCycles, ")! Use parameter 'allowForLarger' to ignore and put read in truncated read group.");
		} else if(alignment.length() >= settings.maxCycles){
			//add to truncated read group
			alignment.setReadGroup(settings.altReadGroupId);
		}

		//add as ready for writing
		addToContainer(_alignmentStorage, &alignment, true);

	} else if(settings.type == ReadGroupType::paired){
		//is orphan
		if(_keepOrphans){
			//add as ready for writing
			addToContainer(_alignmentStorage, &alignment, true);
		} else {
			//filter out (ignore) but write reason to bam log
			_bamFile.filterOut(alignment.name(), alignment.isReverseStrand(), alignment.readGroupId(), alignment.refID());
		}
	}
};

bool TAlignmentSplitMerger::_alignmentCanBeWrittenUnchanged(){
	return  !_recalibrate &&
			!_bamFile.curIsPaired() &&
			_alignmentStorage.empty() &&
			_rgSettings.getType(_bamFile.curReadGroupID())==ReadGroupType::unchanged;
}

//-----------------------------------------
// TOverlapQuantifier
//-----------------------------------------

void TOverlapQuantifier::run(){
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
					UERROR("Mates '", alignment->name(), "' are in different read groups!");
				}

				//calculate overlap and fragment length and add to storage
				size_t overlap = _merger.determineOverlapLength(*alignment, mate->alignment());
				size_t fragmentLength = alignment->fragmentLength();

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
