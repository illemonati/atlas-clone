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

#include "TBitSet.h"
#include "TFile.h"
#include "TGenome.h"
#include "TGenotypeData.h"
#include "TTask.h"

namespace GenomeTasks {

//---------------------------------
// TPileup
//---------------------------------
class TPileup : public TGenome_windows {
private:
	coretools::TOutputFile _out;

	// what to print?
	coretools::TBitSet<8> _printSettings;
	enum {OnlySitesWithData, Depth, Bases, Qualities, Alleles, Mates, Strand, Likelihoods};

	void _handleWindow() override;
	void _handleAlignment() override {}
public:
	TPileup();
	void printPileup();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_pileup : public coretools::TTask {
public:
	TTask_pileup() { _explanation = "Printing pileup from BAM file"; };

	void run() {
		TPileup pileup;
		pileup.printPileup();
	};
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TPILEUP_H_ */
