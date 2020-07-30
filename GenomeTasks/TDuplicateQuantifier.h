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
	TCountDistributionVector _countsPerReadGroup;
	TCountDistribution _countsCombined;

	std::vector<uint32_t> _countsAtPos;
	BAM::TGenomePosition _curPos;
	BAM::TGenomePosition _curChrEnd;

	void _addCurCounts(const BAM::TGenomePosition & nextPos);
	void _handleAlignments();

public:
	TDuplicateQuantifier(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void estimateDuplicationCounts();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_duplicationQuantifier:public TTask{
public:
	TTask_duplicationQuantifier(){ _explanation = "Quantifying read duplication"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TDuplicateQuantifier duplicationQuantifier(Parameters, Logfile, _randomGenerator);
		duplicationQuantifier.estimateDuplicationCounts();
	};
};


}; // end namespace

#endif /* TDUPLICATIONQUANTIFYER_H_ */
