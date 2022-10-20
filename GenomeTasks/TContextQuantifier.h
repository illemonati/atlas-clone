/*
 * TContextQuantifer.h
 *
 *  Created on: Jun 7, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TCONTEXTQUANTIFIER_H_
#define GENOMETASKS_TCONTEXTQUANTIFIER_H_

#include <string>

#include "TGenome.h"
#include "coretools/Main/TTask.h"
#include "coretools/Math/counters.h"

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

#endif /* GENOMETASKS_TCONTEXTQUANTIFIER_H_ */
