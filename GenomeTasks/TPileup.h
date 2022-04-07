/*
 * TPileup.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TPILEUP_H_
#define GENOMETASKS_TPILEUP_H_

#include <array>
#include <set>
#include <string>

#include "TFile.h"
#include "TGenome.h"
#include "TGenotypeData.h"
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TTask.h"

namespace GenomeTasks{


//---------------------------------
// TPileup
//---------------------------------
class TPileup:public TGenome_windows{
private:
	coretools::TOutputFile out;
	bool printOnlySitesWithData;

	//tmp variables
	GenotypeLikelihoods::TBaseCounts _alleleCounts;
	std::array<int, 2> _counts;

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
	TPileup(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void printPileup();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_pileup:public coretools::TTask{
public:
	TTask_pileup(){ _explanation = "Printing pileup from BAM file"; };

	void run(){
		using namespace coretools::instances;
		TPileup pileup(parameters(), &logfile(), &randomGenerator());
		pileup.printPileup();
	};
};


}; // end namespace

#endif /* GENOMETASKS_TPILEUP_H_ */
