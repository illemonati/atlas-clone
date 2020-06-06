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

namespace GenomeTasks{

//----------------------------------------------
// TDuplicateQuantifyer
//----------------------------------------------
class TDuplicateQuantifyer:public TGenome_filtered{
private:
	TCountDistributionVector _countsPerReadGroup;
	TCountDistribution _countsCombined;

	std::vector<uint32_t> _countsAtPos;
	uint32_t _curPos;
	uint32_t _curChrLength;

	void _addCurCounts(const uint32_t nextPos);
	void _handleAlignments();

public:
	TDuplicateQuantifyer(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void estimateDuplicationCounts();
};

}; // end namespace

#endif /* TDUPLICATIONQUANTIFYER_H_ */
