#ifndef TBAMWINDOW_H_
#define TBAMWINDOW_H_

#include "TAlignment.h"
#include "coretools/Main/TError.h"
#include "coretools/Math/TNumericRange.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/TAlleles.h"
#include "genometools/TBed.h"

#include "TSite.h"
#include "genometools/TFastaReader.h"
#include <cstdint>

namespace BAM {

class TBamWindow {
public:
	static constexpr size_t maxReadID = uint16_t(-1);

private:
	genometools::TGenomePosition _from = {size_t(-1), size_t(-1)};
	std::vector<GenotypeLikelihoods::TSite> _entries;
	std::vector<std::vector<uint16_t>> _readIDs;
	std::vector<BAM::TAlignment> _overlap;
	std::vector<bool> _masked;

	genometools::TBed _regions;
	genometools::TAlleles _alleles;

	size_t _iSite      = 0;
	size_t _numMasked  = 0;
	size_t _numBases   = 0;
	size_t _numReads   = 0;
	size_t _sitesData  = 0;
	size_t _sites2Plus = 0;

	coretools::TNumericRange<size_t> _depthFilter{0, true, -1, true};
	coretools::Probability _downProb{1.};
	size_t _upToDepth = 1000;
	bool _filterCpG   = false;
	bool _skipEmpty   = false;
	bool _filterRefN  = false;
	bool _tagReads    = false;

	void _makeBedOrAlleles(const genometools::TChromosomes &Chromosomes, genometools::Morphic Mo);

public:
	TBamWindow(const genometools::TChromosomes &Chromosomes, genometools::Morphic Mo);
	TBamWindow(const genometools::TBed& Regions) : _regions(Regions) {}
	TBamWindow(const genometools::TAlleles& Alleles) : _alleles(Alleles) {}

	genometools::TGenomePosition from() const noexcept {return _from;}
	genometools::TGenomePosition to() const noexcept {return _from + _entries.size();}
	genometools::TGenomePosition position() const noexcept { return from() + _iSite; }
	genometools::TGenomeWindow window() const noexcept {return {from(), to()};}

	size_t size() const noexcept {return _entries.size();}
	bool empty() const noexcept {return _entries.empty();}

	void move(genometools::TGenomeWindow Window, const genometools::TFastaReader& Reference);
	void add(const TAlignment& Alignment);
	void filter();

	void mask(size_t I);
	bool masked(size_t I) const noexcept {return _masked[I];}

	const GenotypeLikelihoods::TSite &operator[](size_t I) const noexcept(coretools::noDebug);
	GenotypeLikelihoods::TSite &operator[](size_t I) noexcept(coretools::noDebug);

	coretools::TConstView<uint16_t> readIDs(size_t I) const noexcept(coretools::noDebug) {
		DEBUG_ASSERT(_tagReads && I < _readIDs.size());
		return _readIDs[I];
	}
	coretools::TConstView<uint16_t> readIDs() const noexcept(coretools::noDebug) { return readIDs(_iSite); }

	const GenotypeLikelihoods::TSite& site() const noexcept(coretools::noDebug) {DEBUG_ASSERT(_iSite < _entries.size()); return _entries[_iSite];}
	GenotypeLikelihoods::TSite& site() noexcept(coretools::noDebug) {DEBUG_ASSERT(_iSite < _entries.size()); return _entries[_iSite];}
	bool endOfSites() const noexcept {return _iSite >= _entries.size();}
	void nextSite();

	size_t numSites() const noexcept { return size() - _numMasked; }
	size_t numReads() const noexcept { return _numReads; }
	size_t numBases() const noexcept { return _numBases; }
	size_t numSitesWithData() const noexcept { return _sitesData; }
	double fracMissing() const noexcept { return (numSites() - numSitesWithData()) / double(numSites()); }
	double depth() const noexcept { return _numBases / double(numSites()); }
	size_t upToDepth() const noexcept { return _upToDepth; }

	const genometools::TBed &regions() const noexcept { return _regions; }
	const genometools::TAlleles &alleles() const noexcept { return _alleles; }

	void requireReadIDs(bool Yes = true) noexcept {_tagReads = Yes;}
	void filterRefN(bool Yes = true) noexcept {_filterRefN = Yes;}
	void skipEmpty(bool Yes=true) noexcept {_skipEmpty = Yes;}
	bool filtersCpG() const noexcept {return _filterCpG;}
};
	
}

#endif
