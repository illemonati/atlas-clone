/*
 * TPMDSCalculator.h
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TPMDSCALCULATOR_H_
#define GENOMETASKS_TPMDSCALCULATOR_H_

#include "TAlignmentTraverser.h"
#include "coretools/Types/probability.h"
#include "coretools/Math/TNumericRange.h"

#include "TOutputBamFile.h"

namespace GenomeTasks{

//----------------------------------------------
// TPMDSCalculator
//----------------------------------------------
class TPMDSCalculator final {
private:
	TAlignmentTraverser _alnTraverser;
	coretools::Probability _pi;
	bool _doFilter = false;
	coretools::TNumericRange<double> _filterRange;
	BAM::TOutputBamFile _outBam;

	std::vector<size_t> _dPMDs_pos;
	std::vector<size_t> _dPMDs_neg;

	double _calculatePMDS(BAM::TAlignment& alignment);
	void _traverseAlignments();

public:
	TPMDSCalculator();
	void run();
};

}; // end namespace


#endif /* GENOMETASKS_TPMDSCALCULATOR_H_ */
