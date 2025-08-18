/*
 * TFinches.h
 *
 *  Created on: Dec 6, 2022
 *      Author: Andreas
 */

#ifndef GENOMETASKS_TFINCHES_H_
#define GENOMETASKS_TFINCHES_H_

#include "TBamTraverser.h"
#include "TOutputBamFile.h"
#include "genometools/TAlleles.h"

namespace BAM {
class TReadGroupMap;
}

namespace GenomeTasks {

struct TPosition {
	genometools::TGenomePosition pos;
	genometools::Base major;
	genometools::Base minor;
TPosition(size_t Ref, size_t Pos, genometools::Base Major, genometools::Base Minor) :
	pos(Ref, Pos), major(Major), minor(Minor) {}
};

//----------------------------------------
// TPMDEstimator.h
//----------------------------------------

class TReadExtractor : public TBamReadTraverser<ReadType::Parsed> {
private:
	size_t _nRef   = 0;
	size_t _nAlt   = 0;
	size_t _nOther = 0;
	size_t _nBoth  = 0;

	BAM::TOutputBamFile _outBAM;
	genometools::TAlleles _alleles;
	void _handleAlignment(BAM::TAlignment& alignment) override;

public:
	TReadExtractor();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TPMDESTIMATOR_H_ */
