/*
 * TReadGroupMerger.h
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TREADGROUPMERGER_H_
#define GENOMETASKS_TREADGROUPMERGER_H_

#include "TGenome.h"

namespace GenomeTasks{

//--------------------------------------
// TReadGroupMerger
//--------------------------------------
class TReadGroupMerger:public TGenome_basic{
private:
	std::vector<uint16_t> readGroupMap;

public:
	TReadGroupMerger(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void mergeReadGroups();
};

}; // end namespace

#endif /* GENOMETASKS_TREADGROUPMERGER_H_ */
