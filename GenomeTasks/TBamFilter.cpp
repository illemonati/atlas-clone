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

void TAlignmentMergerEntry::setAsNonProperPair() {
	_alignment->setIsProperPair(false);
	_ready = true;
};

bool TAlignmentMergerEntry::operator<(const TAlignmentMergerEntry & other) const{
	return _alignment->position() < other._alignment->position();
}

//-----------------------------------------
// TBamFilter
//-----------------------------------------

void TBamFilter::_handleMates(BAM::TAlignment & alignment, iterator mate){
	if(!alignment.isProperPair()){
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
	return	!_recalibrate && 
			!_genome.bamFile().curIsPaired() && 
			_alignmentStorage.empty() &&
			(_removeSoftClippedBases ? (_genome.bamFile().curCIGAR().lengthSoftClippedRight() < _maxNumberOfSoftClippedBases && _genome.bamFile().curCIGAR().lengthSoftClippedLeft() < _maxNumberOfSoftClippedBases) : true);
}

}; //end namespace BamFilter

}; //end namespace





