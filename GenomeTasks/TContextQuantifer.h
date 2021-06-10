/*
 * TContextQuantifer.h
 *
 *  Created on: Jun 7, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TCONTEXTQUANTIFER_H_
#define GENOMETASKS_TCONTEXTQUANTIFER_H_

#include "TGenome.h"
#include "counters.h"
#include "TTask.h"

namespace GenomeTasks{

//----------------------------------------------
// TContextQuantifier
//----------------------------------------------
class TContextQuantifier:public TGenome_parsed{
private:
	coretools::TCountDistributionVector _contextCounts;

	void _handleAlignment();
public:
	TContextQuantifier(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void quantifyContexts();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_quantifyContext:public coretools::TTask{
public:
	TTask_quantifyContext(){ _explanation = "Writing context statistics to file"; };

	void run(coretools::TParameters & Parameters, coretools::TLog* Logfile){
		TContextQuantifier quantifier(Parameters, Logfile, _randomGenerator);
		quantifier.quantifyContexts();
	};
};

}; //end namespace

#endif /* GENOMETASKS_TCONTEXTQUANTIFER_H_ */
