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
#include "coretools/Main/TTask.h"

namespace GenomeTasks{

//--------------------------------------
// TReadGroupMerger
//--------------------------------------
class TReadGroupMerger:public TGenome_basic{
private:
	std::vector<uint16_t> readGroupMap;
public:
	TReadGroupMerger();
	void run();
};

} // end namespace

#endif /* GENOMETASKS_TREADGROUPMERGER_H_ */
