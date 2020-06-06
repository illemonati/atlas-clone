/*
 * TPileup.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TPILEUP_H_
#define GENOMETASKS_TPILEUP_H_

#include "TGenome.h"

namespace GenomeTasks{

//TODO: add flags so users can choose what to write

//---------------------------------
// TPileup
//---------------------------------
class TPileup:public TGenome_windows{
private:
	TOutputFile out;
	bool printOnlySitesWithData;

	//tmp variables
	GenotypeLikelihoods::TBaseData _alleleCounts;
	int _counts[2];

	void _handleWindow();


public:
	TPileup(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void printPileup();
};

}; // end namespace

#endif /* GENOMETASKS_TPILEUP_H_ */
