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

namespace GenomeTask{

//----------------------------------------------
// TContextQuantifier
//----------------------------------------------
class TContextQuantifier:public TGenome_parsed{
private:
	TCountDistributionVector _contextCounts;

	void _handleAlignment();
public:
	TContextQuantifier(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void quantifyContexts();
};


}; //end namespace

#endif /* GENOMETASKS_TCONTEXTQUANTIFER_H_ */
