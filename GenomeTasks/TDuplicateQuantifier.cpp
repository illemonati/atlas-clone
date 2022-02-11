/*
 * simpleGenomeTasks.cpp
 *
 *  Created on: Jun 2, 2020
 *      Author: phaentu
 */


#include "TDuplicateQuantifier.h"

namespace GenomeTasks{

//----------------------------------------------
// TDuplicateQuantifyer
//----------------------------------------------
TDuplicateQuantifier::TDuplicateQuantifier(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator):TGenome_filtered(Parameters, Logfile, RandomGenerator){};

void TDuplicateQuantifier::_addCurCounts(const BAM::TGenomePosition & nextPos){
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

void TDuplicateQuantifier::_handleAlignments(){
	//add to counts
	if(_bamFile.chrChanged()){
		if(_curChrEnd.position() > 0){
			_addCurCounts(_curChrEnd);
		}
		_curPos = _bamFile.curChromosome().chrStart;
		_curChrEnd = _bamFile.curChromosome().chrEnd;
		std::fill(_countsAtPos.begin(), _countsAtPos.end(), 0);
	}

	if(_bamFile.curPosition() == _curPos){
		++_countsAtPos[_bamFile.curReadGroupID()];
	} else if(_bamFile.curPosition() > _curPos){
		_addCurCounts(_bamFile.curPosition());

		//set counts at current position
		_curPos = _bamFile.curPosition();
		std::fill(_countsAtPos.begin(), _countsAtPos.end(), 0);
		_countsAtPos[_bamFile.curReadGroupID()] = 1;
	}
};

void TDuplicateQuantifier::estimateDuplicationCounts(){
	//assembles distribution of how often a read is duplicated
	//now: just how many reads start at the same positions
	_curChrEnd.clear();
	_countsPerReadGroup.resize(_bamFile.numReadGroups());
	_countsAtPos.resize(_bamFile.numReadGroups());

	//iterate through BAM file
	_traverseBAMPassedQC();

	//write output
	std::string filename = _outputName + "_readStartsPerSite.txt";
	_logfile->listFlush("Writing distribution of read starts per site to '" + filename + "' ...");
	coretools::TOutputFile out(filename, {"readGroup", "numReadStarts", "counts"});
	_countsCombined.write(out, "allReadGroups");

	std::vector<std::string> readGroupNames;
	_bamFile.readGroups().fillVectorWithNames(readGroupNames);
	_countsPerReadGroup.write(out, readGroupNames);
	_logfile->done();
};

}; // end namespace
