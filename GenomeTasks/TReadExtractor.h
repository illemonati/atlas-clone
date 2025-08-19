#ifndef GENOMETASKS_TFINCHES_H_
#define GENOMETASKS_TFINCHES_H_

#include "TBamTraverser.h"
#include "TOutputBamFile.h"
#include "genometools/TAlleles.h"

namespace GenomeTasks {

class TReadExtractor : public TBamReadTraverser<ReadType::Parsed> {
private:
	size_t _nRef   = 0;
	size_t _nAlt   = 0;
	size_t _nOther = 0;
	size_t _nBoth  = 0;

	BAM::TOutputBamFile _outBAM;
	genometools::TAlleles _alleles;
	genometools::TAlleles::const_iterator _alIt;
	void _handleAlignment(BAM::TAlignment& alignment) override;

public:
	TReadExtractor();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TPMDESTIMATOR_H_ */
