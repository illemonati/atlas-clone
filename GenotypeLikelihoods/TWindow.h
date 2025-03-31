/*
 * TWindow.h
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#ifndef TWINDOW_H_
#define TWINDOW_H_

#include <cstdint>
#include <string>
#include <vector>

#include "coretools/Types/probability.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/GenomePositions/TGenomeWindow.h"

#include "genometools/Genotypes/Containers.h"
#include "TAlignment.h"
#include "TSite.h"
#include "genometools/TAlleles.h"


namespace coretools {
template<typename T> class TNumericRange;
}
namespace genometools {
class TBed;
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
	size_t _lastReadID = 0;
	mutable std::vector<std::vector<uint32_t>> _readIDs;
	std::vector<BAM::TAlignment> _overlap;
	std::vector<TSite> _sites;
	std::vector<bool> _masked;
	std::string _chrName;

	mutable double _depth                   = 0;
	mutable double _fractionMissing         = 0.;
	mutable double _fractionDepthAtLeastTwo = 0.;
	mutable double _fractionRefIsN          = 0.;
	mutable size_t _numSitesWithData        = 0;

	size_t _numMaskedSites   = 0;
	bool _passedFilters      = false;
	bool _referenceBaseAdded = false;

	void _calcDepth() const;

	// fill sites and clean
	size_t _findFirstPositionWithinWindow(const BAM::TAlignment &Alignment) const;
	void _fillSites(const BAM::TAlignment &alignment);
	void _fillSites(BAM::TAlignment &alignment, const genometools::TAlleles &alleles);
	int _fillSitesDownsampling(std::vector<TSite> &sites, const coretools::Probability &downsamplingProb) const;

	// void stealFromOther(TWindow & other);
	void _downsampleFrom(const TWindow &other, const coretools::Probability &downsamplingProb);

public:
	TWindow(size_t refID, std::string_view ChrName) : genometools::TGenomeWindow(refID, 0), _chrName(ChrName) {}
	TWindow(const TWindow &other, const coretools::Probability &downsamplingProb, size_t UpToDepth, bool Shuffle); 

	// Allow to set chromosome name when jumping
	void move(const genometools::TGenomeWindow &Window);
	void move(const TWindow &Window, std::string_view ChrName);

	void downsampleSites(size_t UpToDepth, bool Shuffle);
	void downsampleSites(coretools::Probability p);

	void addReferenceBaseToSites(const genometools::TAlleles &Alleles);
	void addReferenceBaseToSites(const genometools::TFastaReader &reference);

	size_t applyMask(genometools::TBed &mask, bool doInverseMasking);
	void maskCpG(const genometools::TFastaReader &reference);
	genometools::TBaseProbabilities estimateBaseFrequencies() const;
	void applyDepthFilter(const coretools::TNumericRange<size_t> &DepthRange);

	// getters
	TSite &operator[](size_t internalPos) noexcept { return _sites[internalPos]; }
	const TSite &operator[](size_t internalPos) const noexcept { return _sites[internalPos]; }
	const std::string &chrName() const noexcept { return _chrName; }
	genometools::TGenomePosition position(size_t internalPos) const noexcept { return from() + internalPos; }
	size_t positionOnChr(size_t internalPos) const noexcept { return from().position() + internalPos; }
	bool isMasked(size_t internalPos) const noexcept { return _masked[internalPos]; }

	size_t numReadsInWindow() const noexcept { return _lastReadID; }
	double depth() const noexcept;

	size_t numMaskedSites() const noexcept {
		_calcDepth();
		return _numMaskedSites;
	}
	size_t numSites() const noexcept {return size() - numMaskedSites();}

	size_t numSitesWithData() const noexcept {
		_calcDepth();
		return _numSitesWithData;
	}
	double fracMissing() const noexcept {
		_calcDepth();
		return _fractionMissing;
	}
	void dataSummary() const noexcept;
	bool filter(double maxFracMissing, double maxRefN);
	bool passedFilters() const noexcept { return _passedFilters; }

	// loop over sites
	const std::vector<TSite>& sites() const noexcept {return _sites;}
	iterator begin() noexcept { return _sites.begin(); };
	const_iterator begin() const noexcept { return _sites.begin(); }
	const_iterator cbegin() const noexcept { return _sites.cbegin(); }

	iterator end() noexcept { return _sites.end(); }
	const_iterator end() const noexcept { return _sites.end(); }
	const_iterator cend() const noexcept { return _sites.cend(); }

	void addAlignment(const BAM::TAlignment& usedAlignment);
	void fillSites(const genometools::TAlleles& alleles);
};

} // namespace GenotypeLikelihoods

#endif /* TWINDOW_H_ */
