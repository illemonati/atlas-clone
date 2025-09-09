#ifndef TSITETRAVERSER_H_
#define TSITETRAVERSER_H_

#include "TAlignmentTraverser.h"
#include "TBamWindow.h"
#include "coretools/Containers/TView.h"
#include "coretools/Main/TError.h"
#include "coretools/Types/probability.h"
#include "genometools/TFastaReader.h"
#include <cstdint>

namespace BAM{

class TSiteTraverser {
	std::vector<TAlignmentTraverser> _alnTraversers;

	size_t _wSize = 1'000'000;
	genometools::TGenomeWindow _window;
	genometools::TBed _windowList;
	TBamWindow _bamWindow;

	size_t _i        = 0;
	size_t _iWindows = 0;
	bool _winChanged = true;

	size_t _numBasesTot = 0;
	size_t _numSitesTot = 0;
	size_t _numSDataTot = 0;
	size_t _numReadsTot = 0;
	size_t _numBasesChr = 0;
	size_t _numSitesChr = 0;
	size_t _numSDataChr = 0;
	size_t _numReadsChr = 0;

	coretools::TTimer _timer;
	coretools::TNumericRange<size_t> _depthFilter{0, true, -1, true};
	coretools::Probability _downProb{1.};

	bool _filterCpG   = false;
	bool _filterEmpty = false;
	bool _filterRefN  = false;

	bool _atStart() const noexcept {return _window.size() == 0;}

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
	const genometools::TFastaReader& reference() const noexcept { return _traverser().reference(); };
	const BAM::TBamFile &bamFile(size_t I = 0) const noexcept { return _traverser(I).bamFile(); }
	BAM::TBamFile &bamFile(size_t I = 0) noexcept { return _traverser(I).bamFile(); }

	bool endOfCurChr() const { return _window.from() >= curChr().to(); }
	bool endOfChrs();
	void nextChr();

	size_t numSites() const noexcept {return _numSitesTot;}

	// Per Site access
	void nextSite();
	const GenotypeLikelihoods::TSite& site() const { return _bamWindow[_i]; }
	coretools::TConstView<uint16_t> readIDs() const {return _bamWindow.readIDs(_i);}
	genometools::TGenomePosition position() const noexcept { return _bamWindow.from() + _i; }
	bool winChanged() const noexcept {return _winChanged;}

	// Per Window access
	void nextWindow();
	const TBamWindow& window() const noexcept {return _bamWindow;}

	void requireReference() const;
	void requireSingleBAM() const;
	void requireReadIDs() noexcept {DEV_ASSERT(_atStart()); _bamWindow.requireReadIDs();}

	void filterRefN(bool Yes=true) noexcept {DEV_ASSERT(_atStart()); _filterRefN = Yes;}
	void filterEmpty(bool Yes=true) noexcept {DEV_ASSERT(_atStart()); _filterEmpty = Yes;}

	size_t refID() const noexcept {return _window.refID();}
	const genometools::TChromosomes &chromosomes() const noexcept(coretools::noDebug) { return _traverser().chromosomes(); }
	const genometools::TChromosome &curChr() const noexcept(coretools::noDebug) { return chromosomes()[refID()]; }

	const genometools::TBed &regions() const noexcept { return _bamWindow.regions(); }
	const genometools::TAlleles &alleles() const noexcept { return _bamWindow.alleles(); }

	size_t numSamples() const noexcept { return _alnTraversers.size(); }
	size_t upToDepth() const noexcept { return _bamWindow.upToDepth(); }
};
	
}

#endif
