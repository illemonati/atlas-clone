/*
 * simpleGenomeTasks.h
 *
 *  Created on: Jun 2, 2020
 *      Author: phaentu
 */

#ifndef TDUPLICATIONQUANTIFYER_H_
#define TDUPLICATIONQUANTIFYER_H_


#include "TGenome.h"
#include "counters.h"
#include "TTask.h"

namespace GenomeTasks{

//----------------------------------------------
// TDuplicateQuantifyer
//----------------------------------------------
class TDuplicateQuantifier:public TGenome_filtered{
private:
	coretools::TCountDistributionVector _countsPerReadGroup;
	coretools::TCountDistribution _countsCombined;

	std::vector<uint32_t> _countsAtPos;
	BAM::TGenomePosition _curPos;
	BAM::TGenomePosition _curChrEnd;

	void _addCurCounts(const BAM::TGenomePosition & nextPos);
	void _handleAlignments();

public:
	TDuplicateQuantifier(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void estimateDuplicationCounts();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_duplicationQuantifier:public coretools::TTask{
public:
	TTask_duplicationQuantifier(){ _explanation = "Quantifying read duplication"; };

	void run(){
		using namespace coretools::instances;
		TDuplicateQuantifier duplicationQuantifier(parameters(), &logfile(), &randomGenerator());
		duplicationQuantifier.estimateDuplicationCounts();
	};
};


}; // end namespace

#endif /* TDUPLICATIONQUANTIFYER_H_ */
