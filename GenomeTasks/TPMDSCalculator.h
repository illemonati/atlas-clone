/*
 * TPMDSCalculator.h
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TPMDSCALCULATOR_H_
#define GENOMETASKS_TPMDSCALCULATOR_H_

#include <set>
#include <string>

#include "TBamFile.h"
#include "TOutputBamFile.h"
#include "TGenome_OLD.h"
#include "coretools/Main/TTask.h"
#include "coretools/Types/probability.h"
#include "coretools/Math/TNumericRange.h"

namespace GenomeTasks{

//----------------------------------------------
// TPMDSCalculator
//----------------------------------------------
class TPMDSCalculator:public old::TGenome_parsed{
private:
	coretools::Probability _pi;
	bool _doFilter = false;
	coretools::TNumericRange<double> _filterRange;
	BAM::TOutputBamFile _outBam;

	double _calculatePMDS(BAM::TAlignment& alignment);
	void _handleAlignment(BAM::TAlignment& alignment) override;

public:
	TPMDSCalculator();
	void run();
};

}; // end namespace


#endif /* GENOMETASKS_TPMDSCALCULATOR_H_ */
