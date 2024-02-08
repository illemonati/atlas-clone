/*
 * TContextQuantifer.h
 *
 *  Created on: Jun 7, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TCONTEXTQUANTIFER_H_
#define GENOMETASKS_TCONTEXTQUANTIFER_H_


#include "TBamTraverser.h"

namespace GenomeTasks{

//----------------------------------------------
// TContextQuantifier
//----------------------------------------------
class TContextQuantifier:public TBamReadTraverser<ReadType::Parsed> {
private:
	coretools::TCountDistributionVector<> _contextCounts;

	void _handleAlignment(BAM::TAlignment& alignment) override;
public:
	TContextQuantifier();
	void run();
};

}; //end namespace

#endif /* GENOMETASKS_TCONTEXTQUANTIFER_H_ */
