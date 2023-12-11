/*
 * simpleGenomeTasks.h
 *
 *  Created on: Jun 2, 2020
 *      Author: phaentu
 */

#ifndef TDUPLICATIONQUANTIFYER_H_
#define TDUPLICATIONQUANTIFYER_H_

#include <vector>

#include "coretools/Math/counters.h"

#include "genometools/GenomePositions/TGenomePosition.h"

#include "TBamTraverser.h"

namespace GenomeTasks{

//----------------------------------------------
// TDuplicateQuantifyer
//----------------------------------------------
class TDuplicateQuantifier final : public TBamReadTraverser<ReadType::Filtered> {
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
	void run();
};

}; // end namespace

#endif /* TDUPLICATIONQUANTIFYER_H_ */
