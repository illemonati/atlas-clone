/*
 * TWindow.h
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#ifndef TWINDOW_H_
#define TWINDOW_H_

#include <string>
#include <vector>

#include "TAlignment.h"
#include "GenotypeData.h"
#include "TSite.h"
#include "coretools/Math/TSubsamplePicker.h"
#include "coretools/Types/probability.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/GenomePositions/TGenomeWindow.h"

namespace coretools {
template<typename T> class TNumericRange;
}
namespace genometools {
class TBed;
}
namespace genometools {
class TFastaReader;
}

namespace GenotypeLikelihoods {

//---------------------------------------------------------------
// TWindow
//---------------------------------------------------------------
class TWindow final : public genometools::TGenomeWindow {
public:
	using iterator       = std::vector<TSite>::iterator;
	using const_iterator = std::vector<TSite>::const_iterator;

private:
	// alignment stacks and sites
	std::vector<BAM::TAlignment> _usedAlignments;
	std::vector<TSite> _sites;
	std::string _chrName;

	mutable double _depth                   = 0;
	mutable double _fractionSitesNoData     = 0.;
	mutable double _fractionDepthAtLeastTwo = 0.;
	mutable double _fractionRefIsN          = 0.;
	mutable size_t _numSitesWithData        = 0;

	size_t _numReadsInWindow = 0;
	size_t _numMaskedSites   = 0;
	bool _passedFilters      = false;
	bool _referenceBaseAdded = false;

	void _calcDepth() const;

	// fill sites and clean
	size_t _findFirstPositionWithinWindow(const BAM::TAlignment &alignment) const;
	void _fillSites(const BAM::TAlignment &alignment, std::vector<TSite> &sites, size_t readUpToDepth) const;
	void _fillSites(std::vector<TSite> &sites, size_t readUpToDepth);
	int _fillSitesDownsampling(std::vector<TSite> &sites, size_t readUpToDepth,
	                           const coretools::Probability &downsamplingProb) const;

	//fill sites according to subset: templates to allow for different types of subsets
	template<typename SiteSubsetType>
	void _fillSitesSubset(BAM::TAlignment &alignment, std::vector<TSite> &sites,
	                      coretools::TConstView<SiteSubsetType> thesePos, size_t readUpToDepth) {
		size_t p = _findFirstPositionWithinWindow(alignment);

		// position in window where first one = 0
		// p is at first position of read in window
		auto it = thesePos.begin();
		for (; p < alignment.parsedLength(); ++p) {
			if (alignment.isAlignedAtInternalPos(p) && alignment[p] != genometools::Base::N) {
				size_t internalPos = alignment.positionInRef(p) - from();

				// if read extends past window length
				if (internalPos >= size()) break; // since part of the read maps to next window

				// find position in thesePos
				while (it != thesePos.end() && *it < alignment.positionInRef(p)) ++it;

				if (it != thesePos.end() && *it == alignment.positionInRef(p) &&
				    sites[internalPos].depth() < readUpToDepth) {
					sites[internalPos].add(alignment[p]);
				}
			}
		}
	}

	template<typename SiteSubsetType>
	void _fillSitesSubset(std::vector<TSite> &sites, SiteSubsetType &subset, size_t readUpToDepth){
		sites.resize(size());

		auto thesePositions = subset.getPositionInWindow(*this);

		for(auto & a : _usedAlignments){
			_fillSitesSubset(a, sites, thesePositions, readUpToDepth);
		}
	}

	template<typename SiteSubsetType>
	int _fillSitesSubsetDownsampling(std::vector<TSite> &sites, SiteSubsetType &subset, size_t readUpToDepth,
	                                 const coretools::Probability &downsamplingProb){
		sites.resize(size());

		//get positions that are used
		auto thesePositions = subset.getPositionInWindow(*this);

		int counter = 0;
		for(auto & a : _usedAlignments){
			if(coretools::instances::randomGenerator().getRand() < downsamplingProb){
				_fillSitesSubset(a, sites, thesePositions, readUpToDepth);
				++counter;
			}
		}
		return counter;
	}
	void _clear();

public:
	TWindow(size_t refID, std::string_view ChrName) : genometools::TGenomeWindow(refID, 0), _chrName(ChrName) {};
	TWindow(const TWindow &other, const int readUpToDepth, const coretools::Probability &downsamplingProb) {
		downsampleFromOther(other, readUpToDepth, downsamplingProb);
	}

	// Allow to set chromosome name when jumping
	using genometools::TGenomeWindow::move;
	void move(const genometools::TGenomeWindow &Window);
	void move(const TWindow &Window, std::string_view ChrName);

	// move / expand on same chromosome
	void operator+=(size_t length);
	void resize(size_t newLength);

	// void stealFromOther(TWindow & other);
	void downsampleFromOther(const TWindow &other, size_t readUpToDepth, const coretools::Probability &downsamplingProb);
	template<typename SiteSubsetType>
	void downsampleFromOther(TWindow &other, SiteSubsetType &subset, size_t readUpToDepth,
	                         const coretools::Probability &downsamplingProb){
		//clear and set coordinates
		_clear();
		move(other, other.chrName());

		//fill sites by downsampling, recalculate depth
		_numReadsInWindow = other._fillSitesSubsetDownsampling(_sites, subset, readUpToDepth, downsamplingProb);
		_calcDepth();
	};
	
	template<typename SiteSubsetType>
	void addReferenceBaseToSites(const SiteSubsetType &subset) {
		if (!_referenceBaseAdded && subset.hasPositionsInWindow(*this)) {
			// now only run over sites listed in that window
			const auto thesePositions = subset.getPositionInWindow(*this);
			for (auto &it : thesePositions) {
				size_t pos          = it - from();
				_sites[pos].refBase = it.ref();
			}
			_referenceBaseAdded = true;
		}
	};
	void addReferenceBaseToSites(const genometools::TFastaReader &reference);

	size_t applyMask(genometools::TBed &mask, bool doInverseMasking);
	void maskCpG(const genometools::TFastaReader &reference);
	void downsample(size_t maxDepth, const coretools::TSubsamplePicker &picker);
	GenotypeLikelihoods::TBaseProbabilities estimateBaseFrequencies() const;
	void applyDepthFilter(const coretools::TNumericRange<size_t> &DepthRange);

	// getters
	TSite &operator[](size_t internalPos) noexcept { return _sites[internalPos]; };
	const TSite &operator[](size_t internalPos) const noexcept { return _sites[internalPos]; };
	const std::string &chrName() const noexcept { return _chrName; };
	genometools::TGenomePosition position(size_t internalPos) const noexcept { return from() + internalPos; };
	size_t positionOnChr(size_t internalPos) const noexcept { return from().position() + internalPos; };

	size_t numReadsInWindow() const noexcept { return _numReadsInWindow; };
	double depth() const noexcept;
	size_t numMaskedSites() const noexcept {return _numMaskedSites;}
	void dataSummary() const noexcept;
	bool filter(double maxFracMissing, double maxRefN);
	bool passedFilters() const noexcept { return _passedFilters; };

	// loop over sites
	iterator begin() noexcept { return _sites.begin(); };
	const_iterator begin() const noexcept { return _sites.begin(); };
	const_iterator cbegin() const noexcept { return _sites.cbegin(); };

	iterator end() noexcept { return _sites.end(); };
	const_iterator end() const noexcept { return _sites.end(); };
	const_iterator cend() const noexcept { return _sites.cend(); };

	void addAlignment(const BAM::TAlignment& usedAlignment) {
		_usedAlignments.push_back(usedAlignment);
	}

	void fillSites(size_t readUpToDepth);
	template<typename SiteSubsetType>
	void fillSitesSubset(SiteSubsetType &subset, size_t readUpToDepth){
		_fillSitesSubset(_sites, subset, readUpToDepth);
		_numReadsInWindow = _usedAlignments.size();
	}
};

};     // namespace GenotypeLikelihoods

#endif /* TWINDOW_H_ */
