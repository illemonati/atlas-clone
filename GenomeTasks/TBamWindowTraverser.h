#ifndef TBAMWINDOWTRAVERSER_H_
#define TBAMWINDOWTRAVERSER_H_

#include <memory>
#include <string>
#include <vector>

#include "coretools/TTimer.h"
#include "genometools/BED/TBed.h"

#include "TGenome.h"
#include "TParser.h"
#include "TSiteSubset.h"
#include "TAlignment.h"
#include "TWindow.h"
#include "genometools/GenomePositions/TGenomeWindow.h"

namespace GenomeTasks {

class TBamWindowTraverser {
	std::vector<std::vector<genometools::TGenomeWindow>> _windows;

	// window filters
	double _maxMissing;
	double _maxRefN;

	// window limits
	size_t _limitWindows;
	size_t _skipWindows;

	bool _doMasking;
	genometools::TBed _mask;

	bool _applyDepthFilter;
	bool _filterCpG;

	coretools::TNumericRange<size_t> _depthFilter;
	std::unique_ptr<coretools::TSubsamplePicker> _subsamplePicker;

	// contructor functions
	void _setWindowParameters();
	void _setParsingLimits();
	void _setWindowFilters();
	void _setSiteFilters();
	void _setMasks();

	void _fillAlignments(GenotypeLikelihoods::TWindow &window);

protected:
	TGenome _genome;
	TParser _parser;

	// window params
	size_t _windowSize;

	// mask
	bool _considerRegions;

	// sites
	std::unique_ptr<GenotypeLikelihoods::TSiteSubsetPolymorphic> _subsetPolymoprhic;
	std::unique_ptr<GenotypeLikelihoods::TSiteSubsetMonomorphic> _subsetMonomorphic;

	// site filters
	size_t _readUpToDepth;
	size_t _downsampleDepth;

	void _openSiteSubset(const std::string &filename, bool polymoprhic = true);
	void _applyWindowFilters(GenotypeLikelihoods::TWindow &window);

	void _traverseBAMWindows();
	virtual void _handleWindow(GenotypeLikelihoods::TWindow &window) = 0;
	virtual void _handleChromosome(const genometools::TChromosome &Chr)   = 0;

public:
	TBamWindowTraverser();
};

} // namespace GenomeTasks

#endif /* GENOME_H_ */
