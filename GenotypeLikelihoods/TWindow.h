/*
 * TWindow.h
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#ifndef TWINDOW_H_
#define TWINDOW_H_

#include <iosfwd>
#include <set>
#include <stdint.h>
#include <string>
#include <vector>

#include "TAlignment.h"
#include "TGenotypeData.h"
#include "TSite.h"
#include "TSiteSubset.h"
#include "coretools/Containers/TView.h"
#include "coretools/Math/TSubsamplePicker.h"
#include "coretools/Types/probability.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/TFastaReader.h"

namespace coretools {
class TOutputFile;
}
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
class TWindow : public genometools::TGenomeWindow {
protected:
	// alignment stacks and sites
	std::vector<BAM::TAlignment> usedAlignments;
	std::vector<TSite> _sites;
	
	//details
	size_t _numReadsInWindow;
	std::string _chrName;
	
	// depth
	bool _depthCalculated;
	double _depth, _fractionSitesNoData, _fractionDepthAtLeastTwo, _fractionRefIsN;
	size_t _numSitesWithData;
	void _calcDepth();

	bool _passedFilters;
	bool referenceBaseAdded;

	// fill sites and clean
	size_t _findFirstPositionWithinWindow(const BAM::TAlignment &alignment);
	void _fillSites(BAM::TAlignment &alignment, std::vector<TSite> &sites, size_t readUpToDepth);
	void _fillSites(std::vector<TSite> &sites, size_t readUpToDepth);
	int _fillSitesDownsampling(std::vector<TSite> &sites, size_t readUpToDepth,
	                           const coretools::Probability &downsamplingProb);

	//fill sites according to subset: templates to allow for different types of subsets
	template<typename SiteSubsetType>
	void _fillSitesSubset(BAM::TAlignment &alignment, std::vector<TSite> &sites,
	                      coretools::TConstView<SiteSubsetType> thesePos, size_t readUpToDepth) {
		// genomic position of alignment as seen from window perspective
		size_t p = _findFirstPositionWithinWindow(alignment);

		// position in window where first one = 0
		// p is at first position of read in window
		auto it = thesePos.begin();
		for (; p < alignment.parsedLength(); ++p) {
			if (alignment.isAlignedAtInternalPos(p) && alignment[p] != genometools::Base::N) {
				size_t internalPos = alignment.positionInRef(p) - _from;

				// if read extends past window length
				if (internalPos >= size()) break; // since part of the read maps to next window

				// find position in thesePos
				while (it != thesePos.end() && *it < alignment.positionInRef(p)) ++it;

				// add
				if (it != thesePos.end() && *it == alignment.positionInRef(p) &&
				    sites[internalPos].depth() < readUpToDepth) {
					sites[internalPos].add(alignment[p]);
				}
			}
		}
	};

	template<typename SiteSubsetType>
	void _fillSitesSubset(std::vector<TSite> &sites, SiteSubsetType &subset, size_t readUpToDepth){
		sites.resize(size());

		//get positions that are used
		auto thesePositions = subset.getPositionInWindow(*this);

		//add reads in usedAlignments to sites in window
		for(auto & a : usedAlignments){
			//now fill
			_fillSitesSubset(a, sites, thesePositions, readUpToDepth);
		}
	};

	template<typename SiteSubsetType>
	int _fillSitesSubsetDownsampling(std::vector<TSite> &sites, SiteSubsetType &subset, size_t readUpToDepth,
	                                 const coretools::Probability &downsamplingProb){
		sites.resize(size());

		//get positions that are used
		auto thesePositions = subset.getPositionInWindow(*this);

		//add reads in usedAlignments to sites in window
		int counter = 0;
		for(auto & a : usedAlignments){
			//check if alignment is to be used
			if(coretools::instances::randomGenerator().getRand() < downsamplingProb){
				//now fill
				_fillSitesSubset(a, sites, thesePositions, readUpToDepth);
				++counter;
			}
		}
		return counter;
	};
	void _cleanUpUsedAlignments();
	void _clearAllUsedAlignments();

public:
	TWindow();
	TWindow(TWindow &other, const int readUpToDepth, const coretools::Probability &downsamplingProb);

	// Allow to set chromosome name when jumping
	using genometools::TGenomeWindow::move;
	void move(const genometools::TGenomePosition &From, size_t Length, const std::string ChrName);
	void move(const genometools::TGenomePosition &From, const genometools::TGenomePosition &To,
	                  const std::string ChrName);
	void move(const genometools::TGenomeWindow &Window, const std::string ChrName);
	void setChrName(const std::string ChrName);

	// move / expand on same chromosome
	void operator+=(size_t length);
	void operator-=(size_t length);
	void resize(size_t newLength);

	// void stealFromOther(TWindow & other);
	void downsampleFromOther(TWindow &other, size_t readUpToDepth, const coretools::Probability &downsamplingProb);
	template<typename SiteSubsetType>
	void downsampleFromOther(TWindow &other, SiteSubsetType &subset, size_t readUpToDepth,
	                         const coretools::Probability &downsamplingProb){
		//clear and set coordinates
		clear();
		move(other, other.chrName());

		//fill sites by downsampling, recalculate depth
		_numReadsInWindow = other._fillSitesSubsetDownsampling(_sites, subset, readUpToDepth, downsamplingProb);
		_calcDepth();
	};
	void clear();
	
	template<typename SiteSubsetType>
	void addReferenceBaseToSites(const SiteSubsetType &subset) {
		if (!referenceBaseAdded && subset.hasPositionsInWindow(*this)) {
			// now only run over sites listed in that window
			const auto thesePositions = subset.getPositionInWindow(*this);
			for (auto &it : thesePositions) {
				size_t pos          = it - _from;
				_sites[pos].refBase = it.ref();
			}
			referenceBaseAdded = true;
		}
	};
	void addReferenceBaseToSites(const genometools::TFastaReader &reference);

	void applyMask(genometools::TBed &mask, bool doInverseMasking);
	void maskCpG(const genometools::TFastaReader &reference);
	void downsample(size_t maxDepth, const coretools::TSubsamplePicker &picker);
	GenotypeLikelihoods::TBaseProbabilities estimateBaseFrequencies() const;
	void applyDepthFilter(const coretools::TNumericRange<size_t> &DepthRange);

	// getters
	TSite &operator[](size_t internalPos) { return _sites[internalPos]; };
	const TSite &operator[](size_t internalPos) const { return _sites[internalPos]; };
	const std::string &chrName() const { return _chrName; };
	genometools::TGenomePosition position(size_t internalPos) const { return _from + internalPos; };
	size_t positionOnChr(size_t internalPos) const { return _from.position() + internalPos; };

	size_t numReadsInWindow() const { return _numReadsInWindow; };
	double depth();
	double fractionSitesNoData();
	double fractionDepthAtLeastTwo();
	size_t numSitesWithData();
	double fractionRefIsN();
	void dataSummary();
	bool filter(const double maxFracMissing, const double maxRefN);
	bool passedFilters() const { return _passedFilters; };

	explicit operator std::string() const {
		return fmt::format("{}:{}-{}", chrName(), from().position() + 1, to().position());
	}

	// loop over sites
	std::vector<TSite>::iterator begin() { return _sites.begin(); };
	std::vector<TSite>::iterator end() { return _sites.end(); };
	std::vector<TSite>::const_iterator cbegin() const { return _sites.cbegin(); };
	std::vector<TSite>::const_iterator cend() const { return _sites.cend(); };

	//
	void addAlignment(BAM::TAlignment usedAlignment);
	void fillSites(size_t readUpToDepth);
	template<typename SiteSubsetType>
	void fillSitesSubset(SiteSubsetType &subset, size_t readUpToDepth){
		_fillSitesSubset(_sites, subset, readUpToDepth);
		_numReadsInWindow = usedAlignments.size();
	};
};

};     // namespace GenotypeLikelihoods

#endif /* TWINDOW_H_ */
