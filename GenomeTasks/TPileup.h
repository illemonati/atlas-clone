/*
 * TPileup.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TPILEUP_H_
#define GENOMETASKS_TPILEUP_H_

#include "TGenome.h"
#include "TTask.h"

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
	GenotypeLikelihoods::TBaseCounts _alleleCounts;
	int _counts[2];

	//what to print?
	bool _printDepth;
	bool _printBases;
	bool _printQualities;
	bool _printAlleles;
	bool _printMates;
	bool _printStrand;
	bool _printLikelihoods;

	void _parseField(std::set<std::string> & fields, const std::string tag, bool & flag, const std::string explanation);
	void _handleWindow();
public:
	TPileup(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void printPileup();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_pileup:public TTask{
public:
	TTask_pileup(){ _explanation = "Printing pileup from BAM file"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TPileup pileup(Parameters, Logfile, _randomGenerator);
		pileup.printPileup();
	};
};


}; // end namespace

#endif /* GENOMETASKS_TPILEUP_H_ */
