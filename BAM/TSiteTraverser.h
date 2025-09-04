#ifndef TSITETRAVERSER_H_
#define TSITETRAVERSER_H_

#include "TAlignmentTraverser.h"
#include "TBamWindow.h"
#include "coretools/Main/TError.h"

namespace BAM{

class TSiteTraverser {
	std::vector<TAlignmentTraverser> _alnTraversers;

	size_t _wSize = 1'000'000;
	genometools::TGenomeWindow _window;
	TBamWindow _bamWindow;
	genometools::TBed _windowList;

	size_t _i         = 0;
	size_t _iWindows  = 0;

	size_t _numSites  = 0;
	size_t _nextPrint = 1'000'000;
	coretools::TTimer _timer;

	coretools::TNumericRange<size_t> _depthFilter{1, true, -1, true};
	bool _filterCpG = false;

	void _fillWindow();
	void _filterFindI();
	void _makeHaploDiplo();
	void _advanceWindow();
	void _skipShinkFill();
	void _initChr(size_t RefID);

	void _log();
	const TAlignmentTraverser &_traverser(size_t I = 0) const noexcept(coretools::noDebug) {
		DEBUG_ASSERT(I < _alnTraversers.size());
		return _alnTraversers[I];
	}
	TAlignmentTraverser &_traverser(size_t I = 0) noexcept(coretools::noDebug) {
		DEBUG_ASSERT(I < _alnTraversers.size());
		return _alnTraversers[I];
	}

public:
	TSiteTraverser();

	const std::string &outputName(size_t I = 0) const noexcept { return _traverser(I).outputName(); }
	const GenotypeLikelihoods::TErrorModels &errorModels(size_t I = 0) const noexcept { return _traverser(I).errorModels(); };

	bool endOfCurChr() const { return _window.from() >= curChr().to(); }
	bool endOfChrs();
	void nextChr();

	// Per Site access
	void nextSite();
	const GenotypeLikelihoods::TSite& site() { return _bamWindow[_i]; }
	genometools::TGenomePosition position() const noexcept { return _bamWindow.from() + _i; }

	// Per Window access
	void nextWindow();
	const TBamWindow& window() const noexcept {return _bamWindow;}

	void setDepthFilter(size_t Min, size_t Max = -1) noexcept;
	void requireReference() const;
	void requireSingleBAM() const;

	size_t refID() const noexcept {return _window.refID();}
	const genometools::TChromosomes &chromosomes() const noexcept(coretools::noDebug) { return _traverser().chromosomes(); }
	const genometools::TChromosome &curChr() const noexcept(coretools::noDebug) { return chromosomes()[refID()]; }

	const genometools::TBed &regions() const noexcept { return _bamWindow.regions(); }
	const genometools::TAlleles &alleles() const noexcept { return _bamWindow.alleles(); }
	size_t numSamples() const noexcept { return _alnTraversers.size(); }
};
	
}

#endif
