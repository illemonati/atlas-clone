/*
 * simpleGenomeTasks.cpp
 *
 *  Created on: Jun 2, 2020
 *      Author: phaentu
 */



#include "TDuplicateQuantifier.h"

#include <stddef.h>
#include <algorithm>
#include <cstdint>
#include "TBamFile.h"
#include "TBamFilters.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "coretools/Files/TFile.h"
#include "coretools/Main/TLog.h"
#include "TReadGroups.h"

namespace GenomeTasks{

//----------------------------------------------
// TDuplicateQuantifyer
//----------------------------------------------

void TDuplicateQuantifier::_addCurCounts(const genometools::TGenomePosition & nextPos){
	//add current counts and zero for all positions until nextPos
	uint32_t steps = nextPos - _curPos - 1;
	uint32_t sum = 0;
	for(size_t i=0; i<_countsAtPos.size(); ++i){
		sum += _countsAtPos[i];
		_countsPerReadGroup.add(i, _countsAtPos[i]);
		_countsPerReadGroup.add(i, 0, steps);
	}
	_countsCombined.add(sum);
	_countsCombined.add(sum, steps);
};

void TDuplicateQuantifier::_handleAlignment(){
	//add to counts
	if(_genome.bamFile().chrChanged()){
		if(_curChrEnd.position() > 0){
			_addCurCounts(_curChrEnd);
		}
		_curPos = _genome.bamFile().curChromosome().from();
		_curChrEnd = _genome.bamFile().curChromosome().to();
		std::fill(_countsAtPos.begin(), _countsAtPos.end(), 0);
	}

	if(_genome.bamFile().curPosition() == _curPos){
		++_countsAtPos[_genome.bamFile().curReadGroupID()];
	} else if(_genome.bamFile().curPosition() > _curPos){
		_addCurCounts(_genome.bamFile().curPosition());

		//set counts at current position
		_curPos = _genome.bamFile().curPosition();
		std::fill(_countsAtPos.begin(), _countsAtPos.end(), 0);
		_countsAtPos[_genome.bamFile().curReadGroupID()] = 1;
	}
};

void TDuplicateQuantifier::run(){
	using coretools::instances::logfile;
	//assembles distribution of how often a read is duplicated
	//now: just how many reads start at the same positions
	_curChrEnd.clear();
	_countsPerReadGroup.resize(_genome.bamFile().numReadGroups());
	_countsAtPos.resize(_genome.bamFile().numReadGroups());

	//iterate through BAM file
	_traverseBAMPassedQC();

	//write output
	std::string filename = _genome.outputName() + "_readStartsPerSite.txt";
	logfile().listFlush("Writing distribution of read starts per site to '" + filename + "' ...");
	coretools::TOutputFile out(filename, {"readGroup", "numReadStarts", "counts"});
	_countsCombined.write(out, "allReadGroups");

	std::vector<std::string> readGroupNames;
	_genome.bamFile().readGroups().fillVectorWithNames(readGroupNames);
	_countsPerReadGroup.write(out, readGroupNames);
	logfile().done();
};

}; // end namespace
