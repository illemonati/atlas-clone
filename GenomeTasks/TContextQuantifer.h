/*
 * TContextQuantifer.h
 *
 *  Created on: Jun 7, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TCONTEXTQUANTIFER_H_
#define GENOMETASKS_TCONTEXTQUANTIFER_H_

#include <string>

#include "TGenome.h"
#include "TTask.h"
#include "counters.h"

namespace GenomeTasks{

//----------------------------------------------
// TContextQuantifier
//----------------------------------------------
class TContextQuantifier:public TGenome_parsed {
private:
	coretools::TCountDistributionVector<> _contextCounts;

	void _handleAlignment();
public:
	TContextQuantifier();
	void quantifyContexts();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_quantifyContext:public coretools::TTask{
public:
	TTask_quantifyContext(){ _explanation = "Writing context statistics to file"; };

	void run(){
		TContextQuantifier quantifier;
		quantifier.quantifyContexts();
	};
};

}; //end namespace

#endif /* GENOMETASKS_TCONTEXTQUANTIFER_H_ */
