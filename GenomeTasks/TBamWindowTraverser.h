#ifndef TBAMWINDOWTRAVERSER_H_
#define TBAMWINDOWTRAVERSER_H_

#include "TBamWindows.h"
#include "TGenome.h"
#include "TWindow.h"

namespace GenomeTasks {

class TBamWindowTraverser {
	void _fillAlignments(GenotypeLikelihoods::TWindow &window);

protected:
	TGenome _genome;
	TBamWindows _windows;

	void _traverseBAMWindows();

	virtual void _handleWindow(GenotypeLikelihoods::TWindow &window)    = 0;
	virtual void _handleChromosome(const genometools::TChromosome &Chr) = 0;

public:
	TBamWindowTraverser() : _windows(_genome.bamFile().chromosomes()) {
		_genome.bamFile().setFilters(BAM::TBamFilters{true});
	}
};

} // namespace GenomeTasks

#endif /* GENOME_H_ */
