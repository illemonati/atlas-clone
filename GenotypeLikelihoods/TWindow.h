/*
 * TWindow.h
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#ifndef TWINDOW_H_
#define TWINDOW_H_

#include <stdint.h>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "TAlignment.h"
#include "TGenotypeData.h"
#include "TSite.h"
#include "TSiteSubset.h"
#include "coretools/Math/TSubsamplePicker.h"
#include "coretools/Types/probability.h"
#include "coretools/Containers/TView.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/TFastaReader.h"

namespace coretools { class TLog; }
namespace coretools { class TOutputFile; }
namespace coretools { class TRandomGenerator; }
namespace coretools { template <typename T> class TNumericRange; }
namespace genometools { class TBed; }
namespace genometools { class TFastaReader; }

namespace GenotypeLikelihoods{

//forward declaration to enable copy constructor
class TWindow;

//---------------------------------------------------------------
//TWindow_base
//---------------------------------------------------------------
class TWindow_base:public genometools::TGenomeWindow {
protected:
	std::vector<TSite> _sites;
	size_t _numReadsInWindow;
	std::string _chrName;

	bool referenceBaseAdded;
	//depth
	bool _depthCalculated;
	double _depth, _fractionSitesNoData, _fractionDepthAtLeastTwo, _fractionRefIsN;
	size_t _numSitesWithData;
	void _calcDepth();

	bool _passedFilters;

public:
	//size_t startPos;
	//size_t endPos; //end NOT included in window!
	//size_t length;

	TWindow_base();
	TWindow_base(TWindow & other, const int readUpToDepth, const coretools::Probability & downsamplingProb, coretools::TRandomGenerator* randomGenerator);
	virtual ~TWindow_base();

	//Allow to set chromosome name when jumping
	using genometools::TGenomeWindow::move;
	virtual void move(const genometools::TGenomePosition & From, size_t Length, const std::string ChrName);
	virtual void move(const genometools::TGenomePosition & From, const genometools::TGenomePosition & To, const std::string ChrName);
	virtual void move(const genometools::TGenomeWindow & Window, const std::string ChrName);
	void setChrName(const std::string ChrName);

	//move / expand on same chromosome
	virtual void operator+=(size_t length);
	virtual void operator-=(size_t length);
	virtual void resize(size_t newLength);

	//void stealFromOther(TWindow_base & other);
	void downsampleFromOther(TWindow & other, size_t readUpToDepth, const coretools::Probability & downsamplingProb, coretools::TRandomGenerator* randomGenerator);
	void downsampleFromOther(TWindow & other, TSiteSubset & subset, size_t readUpToDepth, const coretools::Probability & downsamplingProb, coretools::TRandomGenerator* randomGenerator);
	void clear();

	void addReferenceBaseToSites(const genometools::TFastaReader & reference);
	void addReferenceBaseToSites(TSiteSubset & subset);
	void applyMask(genometools::TBed & mask, bool doInverseMasking);
	void maskCpG(const genometools::TFastaReader& reference);
	void downsample(size_t maxDepth, const coretools::TSubsamplePicker & picker);
	GenotypeLikelihoods::TBaseProbabilities estimateBaseFrequencies() const;
	void applyDepthFilter(const coretools::TNumericRange<size_t> & DepthRange);

	//getters
	TSite& operator[](size_t internalPos){ return _sites[internalPos]; };
	const TSite& operator[](size_t internalPos) const { return _sites[internalPos]; };
	const std::string& chrName() const{ return _chrName; };
    genometools::TGenomePosition position(size_t internalPos) const{ return _from + internalPos; };
	size_t positionOnChr(size_t internalPos) const{ return _from.position() + internalPos; };

	size_t numReadsInWindow() const { return _numReadsInWindow; };
	double depth();
	double fractionSitesNoData();
	double fractionDepthAtLeastTwo();
	size_t numSitesWithData();
	double fractionRefIsN();
	void dataSummary();
	bool filter(const double maxFracMissing, const double maxRefN, coretools::TLog* Logfile);
	bool passedFilters() const{ return _passedFilters; };

	explicit operator std::string() const { return fmt::format("{}:{}-{}", chrName(), from().position() + 1, to().position()); }

	//loop over sites
	std::vector<TSite>::iterator begin(){ return _sites.begin(); };
	std::vector<TSite>::iterator end(){ return _sites.end(); };
	std::vector<TSite>::const_iterator cbegin() const{ return _sites.cbegin(); };
	std::vector<TSite>::const_iterator cend() const{ return _sites.cend(); };
};

std::ostream& operator<<(std::ostream& os, const TWindow_base & window);

//---------------------------------------------------------------
//TWindow
//---------------------------------------------------------------
class TWindow:public TWindow_base{
friend class TWindow_base;
private:
	//alignment stacks
	std::vector<BAM::TAlignment> usedAlignments;

	void _cleanUpUsedAlignments();
	void _clearAllUsedAlignments();

	//functions to fill sites from alignments
	size_t _findFirstPositionWithinWindow(const BAM::TAlignment & alignment);

	void _fillSites(BAM::TAlignment & alignment, std::vector<TSite> & sites, size_t readUpToDepth);
	void _fillSites(std::vector<TSite> & sites, size_t readUpToDepth);
	int _fillSitesDownsampling(std::vector<TSite> & sites, size_t readUpToDepth, const coretools::Probability & downsamplingProb, coretools::TRandomGenerator* randomGenerator);

	void _fillSitesSubset(BAM::TAlignment & alignmentIt, std::vector<TSite> & sites, coretools::TConstView<TSiteSubsetSite> thesePos, size_t readUpToDepth);
	void _fillSitesSubset(std::vector<TSite> & sites, TSiteSubset & subset, size_t readUpToDepth);
	int _fillSitesSubsetDownsampling(std::vector<TSite> & sites, TSiteSubset & subset, size_t readUpToDepth, const coretools::Probability & downsamplingProb, coretools::TRandomGenerator* randomGenerator);

public:
	TWindow():TWindow_base(){};

	//Overload moving to take care of alignemnts
	//void move(const size_t RefID, const size_t Start, const size_t End, const std::string ChrName);
	void move(const genometools::TGenomePosition & From, size_t Length, const std::string ChrName);
	void move(const genometools::TGenomePosition & From, const genometools::TGenomePosition & To, const std::string ChrName);
	void move(const genometools::TGenomeWindow & Window, const std::string ChrName);
	void operator+=(size_t length);
	void operator-=(size_t length);
	void resize(size_t newLength);

	void fillSites(size_t readUpToDepth);
	void fillSitesSubset(TSiteSubset & subset, size_t readUpToDepth);

	void addAlignment(BAM::TAlignment usedAlignment);
};

}; //end namespace

#endif /* TWINDOW_H_ */
