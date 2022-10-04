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
#include "genometools/GenomePositions/TGenomePosition.h"
#include "TGenotypeData.h"
#include "TSite.h"
#include "TSiteSubset.h"
#include "coretools/Math/TSubsamplePicker.h"
#include "coretools/Types/probability.h"

namespace BAM { class TFastaBuffer; }
namespace coretools { class TLog; }
namespace coretools { class TOutputFile; }
namespace coretools { class TRandomGenerator; }
namespace coretools { template <typename T> class TNumericRange; }
namespace genometools { class TBed; }

namespace GenotypeLikelihoods{

//forward declaration to enable copy constructor
class TWindow;

//---------------------------------------------------------------
//TWindow_base
//---------------------------------------------------------------
class TWindow_base:public genometools::TGenomeWindow{
protected:
	std::vector<TSite> _sites;
	uint32_t _numReadsInWindow;
	std::string _chrName;

	bool referenceBaseAdded;

	//depth
	bool _depthCalculated;
	double _depth, _fractionSitesNoData, _fractionDepthAtLeastTwo, _fractionRefIsN;
	uint32_t _numSitesWithData;
	void _calcDepth();

	bool _passedFilters;

public:
	//unsigned int startPos;
	//unsigned int endPos; //end NOT included in window!
	//unsigned int length;

	TWindow_base();
	TWindow_base(TWindow & other, const int readUpToDepth, const coretools::Probability & downsamplingProb, coretools::TRandomGenerator* randomGenerator);
	virtual ~TWindow_base();

	//Allow to set chromosome name when jumping
	using genometools::TGenomeWindow::move;
	virtual void move(const genometools::TGenomePosition & From, uint32_t Length, const std::string ChrName);
	virtual void move(const genometools::TGenomePosition & From, const genometools::TGenomePosition & To, const std::string ChrName);
	virtual void move(const genometools::TGenomeWindow & Window, const std::string ChrName);
	void setChrName(const std::string ChrName);

	//move / expand on same chromosome
	virtual void operator+=(uint32_t length);
	virtual void operator-=(uint32_t length);
	virtual void resize(uint32_t newLength);

	//void stealFromOther(TWindow_base & other);
	void downsampleFromOther(TWindow & other, uint32_t readUpToDepth, const coretools::Probability & downsamplingProb, coretools::TRandomGenerator* randomGenerator);
	void downsampleFromOther(TWindow & other, TSiteSubset & subset, uint32_t readUpToDepth, const coretools::Probability & downsamplingProb, coretools::TRandomGenerator* randomGenerator);
	void clear();

	void addReferenceBaseToSites(BAM::TFastaBuffer & reference);
	void addReferenceBaseToSites(TSiteSubset & subset);
	void applyMask(genometools::TBed & mask, bool doInverseMasking);
	void maskCpG(BAM::TFastaBuffer & reference);
	void downsample(uint32_t maxDepth, const coretools::TSubsamplePicker & picker);
	GenotypeLikelihoods::TBaseProbabilities estimateBaseFrequencies() const;
	void applyDepthFilter(const coretools::TNumericRange<uint32_t> & DepthRange);

	//getters
	TSite& operator[](uint32_t internalPos){ return _sites[internalPos]; };
	const TSite& operator[](uint32_t internalPos) const { return _sites[internalPos]; };
	const std::string& chrName() const{ return _chrName; };
    genometools::TGenomePosition position(uint32_t internalPos) const{ return _from + internalPos; };
	uint32_t positionOnChr(uint32_t internalPos) const{ return _from.position() + internalPos; };

	uint32_t numReadsInWindow() const { return _numReadsInWindow; };
	double depth();
	double fractionSitesNoData();
	double fractionDepthAtLeastTwo();
	uint32_t numSitesWithData();
	double fractionRefIsN();
	void dataSummary();
	bool filter(const double maxFracMissing, const double maxRefN, coretools::TLog* Logfile);
	bool passedFilters() const{ return _passedFilters; };

	explicit operator std::string() const { return fmt::format("{}:{}-{}", chrName(), from().position() + 1, to().position()); }

	template<typename Iterator>
	Iterator toBuffer(Iterator b) const { return fmt::format_to(b, "{}:{}-{}", chrName(), from().position() + 1, to().position()); }

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
	uint32_t _findFirstPositionWithinWindow(const BAM::TAlignment & alignment);

	void _fillSites(BAM::TAlignment & alignment, std::vector<TSite> & sites, uint32_t readUpToDepth);
	void _fillSites(std::vector<TSite> & sites, uint32_t readUpToDepth);
	int _fillSitesDownsampling(std::vector<TSite> & sites, uint32_t readUpToDepth, const coretools::Probability & downsamplingProb, coretools::TRandomGenerator* randomGenerator);

	void _fillSitesSubset(BAM::TAlignment & alignmentIt, std::vector<TSite> & sites, std::set<TSiteSubsetSite> & thesePos, unsigned int readUpToDepth);
	void _fillSitesSubset(std::vector<TSite> & sites, TSiteSubset & subset, uint32_t readUpToDepth);
	int _fillSitesSubsetDownsampling(std::vector<TSite> & sites, TSiteSubset & subset, uint32_t readUpToDepth, const coretools::Probability & downsamplingProb, coretools::TRandomGenerator* randomGenerator);

public:
	TWindow():TWindow_base(){};

	//Overload moving to take care of alignemnts
	//void move(const uint32_t RefID, const uint32_t Start, const uint32_t End, const std::string ChrName);
	void move(const genometools::TGenomePosition & From, uint32_t Length, const std::string ChrName);
	void move(const genometools::TGenomePosition & From, const genometools::TGenomePosition & To, const std::string ChrName);
	void move(const genometools::TGenomeWindow & Window, const std::string ChrName);
	void operator+=(uint32_t length);
	void operator-=(uint32_t length);
	void resize(uint32_t newLength);

	void fillSites(uint32_t readUpToDepth);
	void fillSitesSubset(TSiteSubset & subset, uint32_t readUpToDepth);

	void addAlignment(BAM::TAlignment usedAlignment);
};

}; //end namespace

#endif /* TWINDOW_H_ */
