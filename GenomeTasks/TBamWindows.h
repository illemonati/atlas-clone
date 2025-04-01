#ifndef TBAMWINDOWS_H_
#define TBAMWINDOWS_H_

#include <vector>

#include "coretools/Types/probability.h"
#include "genometools/TAlleles.h"
#include "genometools/TBed.h"
#include "genometools/GenomePositions/TGenomeWindow.h"

#include "TParser.h"
#include "TWindow.h"

namespace GenomeTasks {

class TBamWindows {
	std::vector<std::vector<genometools::TGenomeWindow>> _windows;
	TParser _parser;

	// window filters
	double _maxMissing;
	double _maxRefN;

	bool _doMasking;
	genometools::TBed _mask;

	bool _applyDepthFilter;
	bool _filterCpG;
	bool _shuffleSites = false;
	bool _allowDownsampling=false;

	coretools::TNumericRange<size_t> _depthFilter;

	// contructor functions
	void _setWindowParameters(const genometools::TChromosomes& chromosomes);
	void _setWindowFilters();
	void _setSiteFilters();
	void _setMasks(const genometools::TChromosomes& chromosomes);

	// window params
	size_t _windowSize;

	// mask
	bool _considerRegions;

	// sites
	genometools::TAlleles _alleles;

	// site filters
	size_t _upToDepth;
	coretools::Probability _downProb;

public:
	TBamWindows(const genometools::TChromosomes& chromosomes);
	void openSiteSubset(const std::string &filename, const genometools::TChromosomes& chromosomes, bool polymoprhic = true);
	void fillSites(GenotypeLikelihoods::TWindow &window);
	void filter(GenotypeLikelihoods::TWindow &window);

	const TParser& parser() const noexcept {return _parser;}
	void requireReference() const;

	const std::vector<genometools::TGenomeWindow>& operator[](size_t refID) const noexcept {return _windows[refID];}
	size_t uptoDepth() const noexcept {return _upToDepth;}
	bool shuffleSites() const noexcept {return _shuffleSites;}
	bool allowDownsampling() const noexcept {return _allowDownsampling;}
	void allowDownsampling(bool Yes) noexcept {_allowDownsampling = Yes;}
	size_t windowSize() const noexcept {return _windowSize;}
	bool considerRegions() const noexcept {return _considerRegions;}
	const genometools::TAlleles& alleles() const noexcept {return _alleles;}
};
} // namespace GenomeTasks

#endif
