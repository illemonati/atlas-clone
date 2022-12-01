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

#include "coretools/Types/strongTypes.h"

namespace GenomeTasks{

namespace BamFilter{

using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using namespace coretools::str;

//-----------------------------------------
// TAlignmentMergerEntry
//-----------------------------------------

TAlignmentMergerEntry::TAlignmentMergerEntry(BAM::TAlignment* Alignment, bool readyForWriting){
	_alignment = Alignment;
	_ready = readyForWriting;
};

TAlignmentMergerEntry::TAlignmentMergerEntry(TAlignmentMergerEntry && other){
	//copy from other
	_alignment = other._alignment;
	_ready = other._ready;

	//set other to default
	other._alignment = nullptr;
	other._ready = false;
};

TAlignmentMergerEntry::~TAlignmentMergerEntry(){
	delete _alignment;
};

TAlignmentMergerEntry& TAlignmentMergerEntry::operator=(TAlignmentMergerEntry && other){
	if(this != &other){
		//free object
		delete _alignment;

		//copy from other
		_alignment = other._alignment;
		_ready = other._ready;

		//set other to default
		other._alignment = nullptr;
		other._ready = false;
	}

	return *this;
};

const std::string& TAlignmentMergerEntry::name() const {
	return _alignment->name();
};

void TAlignmentMergerEntry::setAsNonProperPair() const{
	_alignment->setIsProperPair(false);
	_ready = true;
};

bool TAlignmentMergerEntry::operator<(const TAlignmentMergerEntry & other) const{
	return _alignment->position() < other._alignment->position();
}

//-----------------------------------------
// TAlignmentStorage
//-----------------------------------------
void addToContainer(TAlignmentStorage & Storage, BAM::TAlignment* Alignment, bool readyForWriting){
	Storage.emplace_back(Alignment, readyForWriting);
}

void addToContainer(TAlignmentStorageSorted & Storage, BAM::TAlignment* Alignment, bool readyForWriting){
	Storage.emplace(Alignment, readyForWriting);
}

//-----------------------------------------
// TBamFilter
//-----------------------------------------
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
			alignment->parse(_genotypeLikelihoodCalculator.sequencingErrorModels());
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
	mate->makeReady();
	_alignmentStorage.emplace_back(&alignment, true);
};

void TBamFilter::_handleSingle(BAM::TAlignment & alignment){
	//read is single end: add for writing
	_alignmentStorage.emplace_back(&alignment, true);
};

bool TBamFilter::_alignmentCanBeWrittenUnchanged(){
	return !_recalibrate && !_bamFile.curIsPaired() && _alignmentStorage.empty();
}

}; //end namespace BamFilter

}; //end namespace





