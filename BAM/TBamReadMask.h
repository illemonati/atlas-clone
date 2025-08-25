#ifndef TBAMREADMASK_H_
#define TBAMREADMASK_H_

#include "coretools/Types/probability.h"
#include "genometools/TBed.h"

namespace BAM {

class TBamReadMask{
private:
	bool _doMasking       = false;
	bool _considerRegions = false;
	bool _maskIsPorous    = false;
	coretools::Probability _porousProb{0.0};
	genometools::TBed _mask;

public:
	void setMasks(const genometools::TChromosomes& Chromosomes);
	[[nodiscard]] bool keepSingle(const genometools::TGenomeWindow& aln) const;
	[[nodiscard]] bool keepPaired(const genometools::TGenomeWindow& aln, const genometools::TGenomeWindow& mate) const;
};
}

#endif
