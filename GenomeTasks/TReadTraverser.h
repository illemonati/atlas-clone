#ifndef TREADTRAVERSER_H_
#define TREADTRAVERSER_H_

#include "TBamFile.h"
#include "TGenome.h"

namespace GenomeTasks {
	
class TReadTraverser {
	TGenome _genome{true};
	bool _eor = false;
public:
	const BAM::TBamFile& bamFile() const noexcept {return _genome.bamFile();}
	BAM::TBamFile& bamFile() noexcept {return _genome.bamFile();}
	const std::string& outputName() const noexcept {return _genome.outputName();}

	const BAM::TRead &read() const noexcept {
		return bamFile().curRead();
	}
	void nextRead() {
		_eor = !bamFile().readNextAlignmentThatPassesFilters();
		if (_eor) {
			bamFile().printEndWithSummary(_genome.outputName());
		} else {
			bamFile().printProgress();
		}
	}
	bool endOfReads() {
		if (bamFile().atStart()) {
			bamFile().startProgressReporting();
			nextRead();
		}
		return _eor;
	}


	const genometools::TChromosomes &chromosomes() const noexcept { return _genome.bamFile().chromosomes(); }
	const genometools::TChromosome &curChr() const noexcept {
		return _genome.bamFile().curChromosome();
	}
};
}

#endif
