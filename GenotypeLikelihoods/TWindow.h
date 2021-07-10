/*
 * TWindow.h
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#ifndef TWINDOW_H_
#define TWINDOW_H_

#include "TBedReaderWindows.h"
#include "TBed.h"
#include "TLog.h"
#include "TSiteSubset.h"
#include "TAlignment.h"
#include "TRandomGenerator.h"
#include "TGenotypeLikelihoodCalculator.h"
#include <vector>
#include "TGenomePosition.h"

namespace GenotypeLikelihoods{

//forward declaration to enable copy constructor
class TWindow;

//---------------------------------------------------------------
//TWindow_base
//---------------------------------------------------------------
class TWindow_base:public BAM::TGenomeWindow{
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

	//TODO: make as much as possible private
public:
	//unsigned int startPos;
	//unsigned int endPos; //end NOT included in window!
	//unsigned int length;

	TWindow_base();
	TWindow_base(TWindow & other, const int readUpToDepth, const Probability & downsamplingProb, TRandomGenerator* randomGenerator);
	virtual ~TWindow_base();

	//Allow to set chromosome name when jumping
	virtual void move(const BAM::TGenomePosition & From, const uint32_t & Length, const std::string ChrName);
	virtual void move(const BAM::TGenomePosition & From, const BAM::TGenomePosition & To, const std::string ChrName);
	virtual void move(const BAM::TGenomeWindow & Window, const std::string ChrName);
	void setChrName(const std::string ChrName);

	//move / expand on same chromosome
	virtual void operator+=(const uint32_t & length);
	virtual void operator-=(const uint32_t & length);
	virtual void resize(const uint32_t & newLength);

	//void stealFromOther(TWindow_base & other);
	void downsampleFromOther(TWindow & other, const uint32_t & readUpToDepth, const Probability & downsamplingProb, TRandomGenerator* randomGenerator);
	void downsampleFromOther(TWindow & other, TSiteSubset & subset, const uint32_t & readUpToDepth, const Probability & downsamplingProb, TRandomGenerator* randomGenerator);
	void clear();

	void addReferenceBaseToSites(BAM::TFastaBuffer & reference);
	void addReferenceBaseToSites(TSiteSubset & subset);
	void applyMask(BAM::TBed & mask, bool doInverseMasking);
	void maskCpG(BAM::TFastaBuffer & reference);
	void downsample(const uint32_t & maxDepth, const coretools::TSubsamplePicker & picker);
	GenotypeLikelihoods::TBaseProbabilities estimateBaseFrequencies() const;
	void applyDepthFilter(const coretools::TNumericRange<uint32_t> & DepthRange);

	//getters
	TSite& operator[](uint32_t internalPos){ return _sites[internalPos]; };
	const TSite& operator[](uint32_t internalPos) const { return _sites[internalPos]; };
	const std::string& chrName() const{ return _chrName; };
	BAM::TGenomePosition position(uint32_t internalPos) const{ return _from + internalPos; };
	uint32_t positionOnChr(uint32_t internalPos) const{ return _from.position() + internalPos; };

	uint32_t numReadsInWindow() const { return _numReadsInWindow; };
	double depth();
	double fractionSitesNoData();
	double fractionDepthAtLeastTwo();
	uint32_t numSitesWithData();
	double fractionRefIsN();
	void dataSummary(TLog* Logfile);
	bool filter(const double maxFracMissing, const double maxRefN, TLog* Logfile);
	bool passedFilters() const{ return _passedFilters; };

	//loop over sites
	std::vector<TSite>::iterator begin(){ return _sites.begin(); };
	std::vector<TSite>::iterator end(){ return _sites.end(); };
	std::vector<TSite>::const_iterator cbegin() const{ return _sites.cbegin(); };
	std::vector<TSite>::const_iterator cend() const{ return _sites.cend(); };
};

std::ostream& operator<<(std::ostream& os, const TWindow_base & window);
coretools::TOutputFile& operator<<(coretools::TOutputFile& out, const TWindow_base & window);

//---------------------------------------------------------------
//TWindow
//---------------------------------------------------------------
class TWindow:public TWindow_base{
friend class TWindow_base;
private:
	//alignment stacks
	std::vector<BAM::TAlignment*> usedAlignments;
	std::vector<BAM::TAlignment*> emptyAlignments;

	void _cleanUpUsedAlignments();
	void _clearAllUsedAlignments();

	//functions to fill sites from alignments
	uint32_t _findFirstPositionWithinWindow(const BAM::TAlignment & alignment);

	void _fillSites(BAM::TAlignment & alignment, std::vector<TSite> & sites, const uint32_t & readUpToDepth);
	void _fillSites(std::vector<TSite> & sites, const uint32_t & readUpToDepth);
	int _fillSitesDownsampling(std::vector<TSite> & sites, const uint32_t & readUpToDepth, const Probability & downsamplingProb, TRandomGenerator* randomGenerator);

	void _fillSitesSubset(BAM::TAlignment & alignmentIt, std::vector<TSite> & sites, std::set<TSiteSubsetSite> & thesePos, const unsigned int & readUpToDepth);
	void _fillSitesSubset(std::vector<TSite> & sites, TSiteSubset & subset, const uint32_t & readUpToDepth);
	int _fillSitesSubsetDownsampling(std::vector<TSite> & sites, TSiteSubset & subset, const uint32_t & readUpToDepth, const Probability & downsamplingProb, TRandomGenerator* randomGenerator);

public:
	TWindow();
	~TWindow();

	//Overload moving to take care of alignemnts
	//void move(const uint32_t RefID, const uint32_t Start, const uint32_t End, const std::string ChrName);
	void move(const BAM::TGenomePosition & From, const uint32_t & Length, const std::string ChrName);
	void move(const BAM::TGenomePosition & From, const BAM::TGenomePosition & To, const std::string ChrName);
	void move(const BAM::TGenomeWindow & Window, const std::string ChrName);
	void operator+=(const uint32_t & length);
	void operator-=(const uint32_t & length);
	void resize(const uint32_t & newLength);

	void review();
	void printStacks();

	void fillSites(const uint32_t & readUpToDepth);
	void fillSitesSubset(TSiteSubset & subset, const uint32_t & readUpToDepth);

	BAM::TAlignment* swapUsedForEmptyAlignment(BAM::TAlignment* usedAlignment);
};

}; //end namespace

#endif /* TWINDOW_H_ */
