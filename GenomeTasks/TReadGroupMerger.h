/*
 * TReadGroupMerger.h
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TREADGROUPMERGER_H_
#define GENOMETASKS_TREADGROUPMERGER_H_

#include "TGenome.h"
#include "TTask.h"

namespace GenomeTasks{

//--------------------------------------
// TReadGroupMerger
//--------------------------------------
class TReadGroupMerger:public TGenome_basic{
private:
	std::vector<uint16_t> readGroupMap;

public:
	TReadGroupMerger(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void mergeReadGroups();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_mergeReadGroups:public TTask{
public:
	TTask_mergeReadGroups(){ _explanation = "Merging read groups in a BAM file"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TReadGroupMerger readGroupMerger(Parameters, Logfile, _randomGenerator);
		readGroupMerger.mergeReadGroups();
	};
};

}; // end namespace

#endif /* GENOMETASKS_TREADGROUPMERGER_H_ */
