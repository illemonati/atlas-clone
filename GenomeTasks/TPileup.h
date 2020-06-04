/*
 * TPileup.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TPILEUP_H_
#define GENOMETASKS_TPILEUP_H_

#include "TGenome.h"

//---------------------------------
// TPileup
//---------------------------------
class TPileup:public TGenome_windows{
private:
	TOutputFile out;
	bool printOnlySitesWithData;

	//tmp variables
	GenotypeLikelihoods::TBaseData _alleleCounts;
	GenotypeLikelihoods::TGenotypeLikelihoods _genoLik;
	int _counts[2];

	void _handleWindow();


public:
	TPileup(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void printPileup();
};



#endif /* GENOMETASKS_TPILEUP_H_ */
