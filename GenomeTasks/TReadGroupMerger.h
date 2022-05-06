/*
 * TReadGroupMerger.h
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TREADGROUPMERGER_H_
#define GENOMETASKS_TREADGROUPMERGER_H_

#include <stdint.h>
#include <string>
#include <vector>

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
	TReadGroupMerger();
	void mergeReadGroups();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_mergeReadGroups:public coretools::TTask{
public:
	TTask_mergeReadGroups(){ _explanation = "Merging read groups in a BAM file"; };

	void run(){
		TReadGroupMerger readGroupMerger;
		readGroupMerger.mergeReadGroups();
	};
};

}; // end namespace

#endif /* GENOMETASKS_TREADGROUPMERGER_H_ */
