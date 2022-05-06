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
	void _handleWindow() override;
	void _handleAlignment() override {}
public:
	TPileup();
	void printPileup();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_pileup:public coretools::TTask{
public:
	TTask_pileup(){ _explanation = "Printing pileup from BAM file"; };

	void run(){
		TPileup pileup;
		pileup.printPileup();
	};
};


}; // end namespace

#endif /* GENOMETASKS_TPILEUP_H_ */
