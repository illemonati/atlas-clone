/*
 * TContextQuantifer.h
 *
 *  Created on: Jun 7, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TCONTEXTQUANTIFER_H_
#define GENOMETASKS_TCONTEXTQUANTIFER_H_

#include "coretools/Math/counters.h"

#include "TBamTraverser.h"

namespace GenomeTasks{

//----------------------------------------------
// TContextQuantifier
//----------------------------------------------
class TContextQuantifier:public TBamTraverser<true> {
private:
	coretools::TCountDistributionVector<> _contextCounts;

	void _handleAlignment(BAM::TAlignment& alignment) override;
public:
	TContextQuantifier();
	void run();
};

}; //end namespace

#endif /* GENOMETASKS_TCONTEXTQUANTIFER_H_ */
