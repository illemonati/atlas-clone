/*
 * simpleGenomeTasks.cpp
 *
 *  Created on: Jun 2, 2020
 *      Author: phaentu
 */


#include <TDuplicateQuantifier.h>

namespace GenomeTasks{

//----------------------------------------------
// TDuplicateQuantifyer
//----------------------------------------------
TDuplicateQuantifier::TDuplicateQuantifier(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator){
	_curChrLength = 0;
	_curPos = 0;
};

void TDuplicateQuantifier::_addCurCounts(const uint32_t nextPos){
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
	if(_bamFile.chrChanged()){
		//write last
	}

	//add to counts
	if(_bamFile.chrChanged()){
		if(_curChrLength > 0){
			_addCurCounts(_curChrLength);
		}
		_curPos = 0;
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
	_curChrLength = 0;
	_curPos = 0;
	_countsPerReadGroup.resize(_bamFile._readGroups.size());
	_countsAtPos.resize(_bamFile._readGroups.size());

	//iterate through BAM file
	_traverseBAMPassedQC();

	//write output
	std::string filename = _outputName + "_readStartsPerSite.txt";
	_logfile->listFlush("Writing distribution of read starts per site to '" + filename + "' ...");
	TOutputFile out(filename, {"readGroup", "numReadStarts", "counts"});
	_countsCombined.write(out, "allReadGroups");

	std::vector<std::string> readGroupNames;
	_bamFile._readGroups.fillVectorWithNames(readGroupNames);
	_countsPerReadGroup.write(out, readGroupNames);
	_logfile->done();
};

}; // end namespace
