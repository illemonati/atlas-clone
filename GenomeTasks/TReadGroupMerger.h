/*
 * TReadGroupMerger.h
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TREADGROUPMERGER_H_
#define GENOMETASKS_TREADGROUPMERGER_H_

#include <vector>

#include "TReadTraverser.h"

namespace GenomeTasks{

//--------------------------------------
// TReadGroupMerger
//--------------------------------------
class TReadGroupMerger {
private:
	TReadTraverser _genome{false};
	std::vector<size_t> _readGroupMap;
public:
	TReadGroupMerger();
	void run();
};

} // end namespace

#endif /* GENOMETASKS_TREADGROUPMERGER_H_ */
