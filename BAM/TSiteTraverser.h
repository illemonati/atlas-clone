#ifndef TSITEGLFTRAVERSER_H_
#define TSITEGLFTRAVERSER_H_

#include "TAlignmentTraverser.h"
#include "TBamWindow.h"

namespace BAM{

class TSiteTraverser {
	TAlignmentTraverser _alnTraverser;

	static constexpr size_t _wSize = 1;
	genometools::TGenomeWindow _window{0, 0, _wSize};
	TBamWindow _bamWindow;
	genometools::TBed _windowList;

	size_t _i         = 0;
	size_t _iWindows  = 0;
	size_t _minDepth  = 1;

	void _fillWindow();
	void _findFirstI();
	void _makeHaploDiplo();
	void _advanceWindow();
	void _skipShinkFill();
	void _initChr(size_t RefID);

public:
	TSiteTraverser();

	const std::string &outputName() const noexcept { return _alnTraverser.outputName(); }
	const GenotypeLikelihoods::TErrorModels &errorModels() const noexcept { return _alnTraverser.errorModels(); };

	bool endOfCurChr() const { return _window.from() >= curChr().to(); }
	bool endOfChrs() const { return refID() >= chromosomes().size(); }
	void nextChr();

	// Per Site access
	void nextSite();
	const GenotypeLikelihoods::TSite& site() { return _bamWindow[_i]; }
	genometools::TGenomePosition position() const noexcept { return _bamWindow.from() + _i; }

	// Per Window access
	void nextWindow();
	const TBamWindow& window() const noexcept {return _bamWindow;}

	void setMinDepth(size_t Depth) noexcept;

	size_t refID() const noexcept {return _window.refID();}
	const genometools::TChromosomes &chromosomes() const noexcept(coretools::noDebug) { return _alnTraverser.chromosomes(); }
	const genometools::TChromosome &curChr() const noexcept(coretools::noDebug) { return chromosomes()[refID()]; }

	const genometools::TBed &regions() const noexcept { return _bamWindow.regions(); }
	const genometools::TAlleles &alleles() const noexcept { return _bamWindow.alleles(); }
};
	
}

#endif
