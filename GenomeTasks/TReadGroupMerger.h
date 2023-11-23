/*
 * TReadGroupMerger.h
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TREADGROUPMERGER_H_
#define GENOMETASKS_TREADGROUPMERGER_H_

#include <vector>

#include "TGenome.h"

namespace GenomeTasks{

//--------------------------------------
// TReadGroupMerger
//--------------------------------------
class TReadGroupMerger {
private:
	TGenome _genome;
	std::vector<uint16_t> readGroupMap;
public:
	TReadGroupMerger();
	void run();
};

} // end namespace

#endif /* GENOMETASKS_TREADGROUPMERGER_H_ */
