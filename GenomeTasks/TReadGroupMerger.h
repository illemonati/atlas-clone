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
	TReadGroupMerger(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void mergeReadGroups();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_mergeReadGroups:public coretools::TTask{
public:
	TTask_mergeReadGroups(){ _explanation = "Merging read groups in a BAM file"; };

	void run(coretools::TParameters & Parameters, coretools::TLog* Logfile){
		TReadGroupMerger readGroupMerger(Parameters, Logfile, _randomGenerator);
		readGroupMerger.mergeReadGroups();
	};
};

}; // end namespace

#endif /* GENOMETASKS_TREADGROUPMERGER_H_ */
