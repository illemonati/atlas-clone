#ifndef TBAMWINDOW_H_
#define TBAMWINDOW_H_

#include "TAlignment.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/TAlleles.h"
#include "genometools/TBed.h"

#include "TSite.h"
#include "genometools/TFastaReader.h"

namespace BAM {

class TBamWindow {
private:
	genometools::TGenomePosition _from = {size_t(-1), size_t(-1)};
	std::vector<GenotypeLikelihoods::TSite> _entries;
	std::vector<BAM::TAlignment> _overlap;
	std::vector<bool> _masked;

	genometools::TBed _regions;
	genometools::TAlleles _alleles;

	size_t _numMasked   = 0;
	size_t _depthTot    = 0;
	size_t _sitesData   = 0;
	size_t _sites2Plus = 0;

	size_t _upToDepth = 1000;

	void _makeBedOrAlleles(const genometools::TChromosomes &Chromosomes);

public:
	TBamWindow(const genometools::TChromosomes &Chromosomes);
	TBamWindow(const genometools::TBed& Regions) : _regions(Regions) {}
	TBamWindow(const genometools::TAlleles& Alleles) : _alleles(Alleles) {}

	genometools::TGenomePosition from() const noexcept {return _from;}
	genometools::TGenomePosition to() const noexcept {return _from + _entries.size();}
	genometools::TGenomeWindow window() const noexcept {return {from(), to()};}

	size_t size() const noexcept {return _entries.size();}
	bool empty() const noexcept {return _entries.empty();}

	void move(genometools::TGenomeWindow Window, const genometools::TFastaReader& Reference, bool FilterCpG);
	void add(const TAlignment& Alignment);
	void mask(size_t I);

	const GenotypeLikelihoods::TSite& operator[](size_t I) const noexcept(coretools::noDebug) {
		DEBUG_ASSERT(I < size());
		return _entries[I];
	}

	size_t numSites() const noexcept { return size() - _numMasked; }
	size_t numSitesWithData() const noexcept { return _sitesData; }
	double fracMissing() const noexcept { return (numSites() - numSitesWithData()) / double(numSites()); }
	double depth() const noexcept { return _depthTot / double(numSites()); }

	const genometools::TBed &regions() const noexcept { return _regions; }
	const genometools::TAlleles &alleles() const noexcept { return _alleles; }

	auto begin() noexcept { return _entries.begin(); }
	auto begin() const noexcept {return _entries.begin();}
	auto cbegin() const noexcept {return _entries.cbegin();}
	auto rbegin() noexcept {return _entries.rbegin();}
	auto rbegin() const noexcept {return _entries.rbegin();}
	auto crbegin() const noexcept {return _entries.crbegin();}
	auto end() noexcept {return _entries.end();}
	auto end() const noexcept {return _entries.end();}
	auto cend() const noexcept {return _entries.cend();}
	auto rend() noexcept {return _entries.rend();}
	auto rend() const noexcept {return _entries.rend();}
	auto crend() const noexcept {return _entries.crend();}
};
	
}

#endif
