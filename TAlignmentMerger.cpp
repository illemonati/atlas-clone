/*
 * TAlignmentMerger.cpp
 *
 *  Created on: Mar 28, 2019
 *      Author: wegmannd
 */

#include "TAlignmentMerger.h"

TAlignmentMerger::TAlignmentMerger(BamTools::BamWriter* Writer, TAlignmentParser* Parser, int MaxDistanceBetweenMates){
	writer = Writer;
	parser = Parser;
	_maxDistanceBetweenMates = MaxDistanceBetweenMates;
	_adaptQuality = true;
	_filterOrphans = true;
};

void TAlignmentMerger::_writeAlignment(std::vector< TAlignmentMergerEntry >::iterator & it){
	//save the alignment to the bam file
	it->alignment->save(*writer, parser->genoMap, parser->minQualForPrinting, parser->maxQualForPrinting, parser->qualMap);
	it = alignmentStorage.erase(it);
};

void TAlignmentMerger::_writeAllThatAreReady(){
	std::vector< TAlignmentMergerEntry >::iterator it = alignmentStorage.begin();

	while(it != alignmentStorage.end() && it->ready){
		_writeAlignment(it);
	}
};

std::vector< TAlignmentMergerEntry >::iterator TAlignmentMerger::_findMate(TAlignment & alignment){
	std::vector< TAlignmentMergerEntry >::iterator it;
	for(it=alignmentStorage.begin(); it!=alignmentStorage.end(); ++it){
		//found its mate!
		if(it->alignment->alignmentName == alignment.alignmentName){
			return it;
		}
	}

	return alignmentStorage.end();
};

void TAlignmentMerger::addToBeMerged(TAlignment & alignment){
	std::vector< TAlignmentMergerEntry >::iterator it = _findMate(alignment);

	if(it == alignmentStorage.end()){
		//no mate found: add to storage
		alignmentStorage.emplace_back(alignment, false);
	} else {
		//mate found, merge!
		parser->mergeAlignedBasesBamReads(it->alignment, &alignment, _adaptQuality);
		it->ready = true;
		alignmentStorage.emplace_back(alignment, true);
	}
};

void TAlignmentMerger::addReadyToBeWritten(TAlignment & alignment){
	if(alignmentStorage.empty()){
		alignment.save(*writer, parser->genoMap, parser->minQualForPrinting, parser->maxQualForPrinting, parser->qualMap);
	} else {
		alignmentStorage.emplace_back(alignment, true);
	}
};

void TAlignmentMerger::addAsImproperPair(TAlignment & alignment){
	if(_filterOrphans){
		//no need to keep mate in list anymore
		parser->removeFromBlacklist(alignment, "not a proper pair (orphan)");
	} else {
		//set to improper read
		alignment.isProperPair = false;
		addReadyToBeWritten(alignment);
	}
};

void TAlignmentMerger::writeUpTo(const int position){
	//writes all that are ready or too far away
	std::vector< TAlignmentMergerEntry >::iterator it = alignmentStorage.begin();

	while(it != alignmentStorage.end() && (it->ready || position - it->alignment->position > _maxDistanceBetweenMates)){
		if(it->ready){
			_writeAlignment(it);
		} else {
			if(_filterOrphans){
				parser->addToBlacklist(*(alignmentStorage.begin()->alignment), "orphaned read: mate is farther away than " + toString(_maxDistanceBetweenMates) + " bp");
				it = alignmentStorage.erase(alignmentStorage.begin());
			} else {
				it->setAsNonProperPair();
				_writeAlignment(it);
			}
		}
	}
};

void TAlignmentMerger::clear(){
	//write everything and mark reads with missing mates as improper.
	std::vector< TAlignmentMergerEntry >::iterator it = alignmentStorage.begin();

	//addToBeMerged reads in storage to blacklist
	while(it != alignmentStorage.end()){
		//
		if(!it->ready){
			if(_filterOrphans){
				//addToBeMerged reads in storage that have not found mate to blacklist and clear from storage
				parser->addToBlacklist(*(it->alignment), "orphaned at chromosome switch");
				it = alignmentStorage.erase(it);
			} else {
				//set reads in storage to improper pairs but ready for writing
				it->setAsNonProperPair();
			}
		}

		//save the alignment to the bam file
		_writeAlignment(it);
	}
};




