/*
 * simpleGenomeTasks.h
 *
 *  Created on: Jun 2, 2020
 *      Author: phaentu
 */

#ifndef TDUPLICATIONQUANTIFYER_H_
#define TDUPLICATIONQUANTIFYER_H_

#include <stdint.h>
#include <string>
#include <vector>

#include "TGenome_OLD.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "coretools/Main/TTask.h"
#include "coretools/Math/counters.h"

namespace GenomeTasks{

//----------------------------------------------
// TDuplicateQuantifyer
//----------------------------------------------
class TDuplicateQuantifier:public old::TGenome_filtered{
private:
	TGenome _genome;
	coretools::TCountDistributionVector<> _countsPerReadGroup;
	coretools::TCountDistribution<> _countsCombined;

	std::vector<uint32_t> _countsAtPos;
    genometools::TGenomePosition _curPos;
    genometools::TGenomePosition _curChrEnd;

	void _addCurCounts(const genometools::TGenomePosition & nextPos);
	void _handleAlignment() override;

public:
	void estimateDuplicationCounts();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_duplicationQuantifier:public coretools::TTask{
public:
	TTask_duplicationQuantifier(){ _explanation = "Quantifying read duplication"; };

	void run(){
		TDuplicateQuantifier duplicationQuantifier;
		duplicationQuantifier.estimateDuplicationCounts();
	};
};


}; // end namespace

#endif /* TDUPLICATIONQUANTIFYER_H_ */
