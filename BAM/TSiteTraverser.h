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

	size_t _iWindows = 0;
	bool _winChanged = true;

	size_t _numBasesTot    = 0;
	size_t _numSitesTot    = 0;
	size_t _numWithDataTot = 0;
	size_t _numReadsTot    = 0;
	size_t _numBasesChr    = 0;
	size_t _numSitesChr    = 0;
	size_t _numWithDataChr = 0;
	size_t _numReadsChr    = 0;

	coretools::TTimer _timer;

	bool _atStart() const noexcept {return _window.size() == 0;}

	void _fillWindow();
	void _makeHaploDiplo();
	void _skipShinkFill();
	void _initChr(size_t RefID);
	void _advanceWindow();

	void _logSites() const;
	void _logWindows() const;

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
	const BAM::RGInfo::TReadGroupInfo &rgInfo(size_t I = 0) const noexcept { return _traverser(I).rgInfo(); }
	BAM::RGInfo::TReadGroupInfo &rgInfo(size_t I = 0) noexcept { return _traverser(I).rgInfo(); }
	const genometools::TFastaReader& reference() const noexcept { return _traverser().reference(); };
	const BAM::TBamFile &bamFile(size_t I = 0) const noexcept { return _traverser(I).bamFile(); }
	BAM::TBamFile &bamFile(size_t I = 0) noexcept { return _traverser(I).bamFile(); }

	bool endOfCurChr() const { return _window.from() >= curChr().to(); }
	bool endOfChrs();
	void nextChr();

	// useful after traversing
	size_t numBases() const noexcept {return _numBasesTot;}
	size_t numSites() const noexcept {return _numSitesTot;}
	size_t numSitesWithData() const noexcept {return _numWithDataTot;}
	size_t numReads() const noexcept {return _numReadsTot;}

	// Per Site access
	void nextSite();
	const GenotypeLikelihoods::TSite& site() const { return _bamWindow.site(); }
	coretools::TConstView<uint16_t> readIDs() const {return _bamWindow.readIDs();}
	genometools::TGenomePosition position() const noexcept { return _bamWindow.position(); }
	bool winChanged() const noexcept {return _winChanged;}

	// Per Window access
	void nextWindow();
	TBamWindow& window() noexcept;

	void requireReference() const;
	void requireSingleBAM() const;
	void requireReadIDs() noexcept {DEV_ASSERT(_atStart()); _bamWindow.requireReadIDs();}

	void filterRefN(bool Yes=true) noexcept;
	void skipEmpty(bool Yes=true) noexcept {DEV_ASSERT(_atStart()); _bamWindow.skipEmpty(Yes);}

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
