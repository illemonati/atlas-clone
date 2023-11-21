#ifndef TBAMWINDOWTRAVERSER_H_
#define TBAMWINDOWTRAVERSER_H_

#include <memory>
#include <string>
#include <vector>

#include "coretools/TTimer.h"
#include "genometools/BED/TBed.h"
#include "genometools/GenomePositions/TChromosomes.h"

#include "TGenome.h"
#include "TParser.h"

#include "TAlignment.h"
#include "TWindow.h"

namespace GenomeTasks {

class TBamWindowTraverser {
	BAM::TAlignment _curAlignment;
	bool _hasWindowIndent;
	coretools::TTimer _windowTimer;

	// contructor functions
	void _setWindowParameters();
	void _setParsingLimits();
	void _setWindowFilters();
	void _setSiteFilters();
	void _setMasks();

	void _jumpToEnd();
	void _setCountersBeginningOfChromosome();
	bool _incrementWindow(GenotypeLikelihoods::TWindow &window);
	bool _moveToNextWindow(GenotypeLikelihoods::TWindow &window);
	bool _incrementPredefinedWindow();
	bool _moveToNextPredefinedWindow(GenotypeLikelihoods::TWindow &window);

	bool _moveWindow(GenotypeLikelihoods::TWindow &window);
	void _readAlignmentsIntoWindow(GenotypeLikelihoods::TWindow &window);
	bool _readAndParseAlignment();
	bool _readDataInNextWindow(GenotypeLikelihoods::TWindow &window);

protected:
	TGenome _genome;
	TParser _parser;
	const genometools::TChromosomes &_chromosomes;
	std::vector<genometools::TChromosome>::const_iterator _curChromosome;

	// window params
	size_t _windowSize;
	size_t _numWindowsOnChr;
	size_t _windowNumber;
	bool _chrChangedWindow;

	// predefined windows
	genometools::TGenomeWindowList _predefinedWindows;
	std::multiset<genometools::TGenomeWindow>::iterator _curPredefinedWindow;

	// window limits
	size_t _limitWindows;
	size_t _skipWindows;

	// window filters
	double _maxMissing;
	double _maxRefN;

	// mask
	bool _doMasking, _considerRegions;
	genometools::TBed _mask;

	// sites
	std::unique_ptr<GenotypeLikelihoods::TSiteSubsetPolymorphic> _subsetPolymoprhic;
	std::unique_ptr<GenotypeLikelihoods::TSiteSubsetMonomorphic> _subsetMonomorphic;

	// site filters
	bool _applyDepthFilter;
	size_t _readUpToDepth;
	coretools::TNumericRange<size_t> _depthFilter;
	bool _filterCpG;
	size_t _downsampleDepth;
	std::unique_ptr<coretools::TSubsamplePicker> subsamplePicker;

	// tmp variables


	void _openSiteSubset(const std::string &filename, bool polymoprhic = true);
	void _applyWindowFilters(GenotypeLikelihoods::TWindow &window);

	void _traverseBAMWindows();
	virtual void _handleWindow(GenotypeLikelihoods::TWindow& window) = 0;

public:
	TBamWindowTraverser();
};

} // namespace GenomeTasks

#endif /* GENOME_H_ */
