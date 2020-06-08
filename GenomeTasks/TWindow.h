/*
 * TWindow.h
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#ifndef TWINDOW_H_
#define TWINDOW_H_

#include "TBedReaderWindows.h"
#include "TLog.h"
#include "TSiteSubset.h"
#include "TAlignment.h"
#include "TQualityMap.h"
#include "TRandomGenerator.h"
#include "TGenotypeLikelihoodCalculator.h"
#include <vector>

namespace GenomeTasks{

//forward declaration to enable copy constructor
class TWindow;

//---------------------------------------------------------------
//TWindow_base
//---------------------------------------------------------------
class TWindow_base{
protected:
	std::vector<GenotypeLikelihoods::TSite> _sites;
	uint32_t _numReadsInWindow;
	uint32_t _refId;
	std::string _chrName;

	TGenotypeMap genoMap;
	bool referenceBaseAdded;

	//depth
	bool _depthCalculated;
	double _depth, _fractionSitesNoData, _fractionDepthAtLeastTwo, _fractionRefIsN;
	uint32_t _numSitesWithData;
	void _calcDepth();

	bool _passedFilters;

	void _setCoordinates(long Start, long End, int RefId);

	//TODO: make as much as possible private
public:
	unsigned int startPos;
	unsigned int endPos; //end NOT included in window!
	unsigned int length;

	TWindow_base();
	TWindow_base(TWindow_base & other);
	virtual ~TWindow_base();

	void stealFromOther(TWindow_base & other);
	TWindow_base(TWindow & other, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);
	void downsampleFromOther(TWindow & other, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);
	void downsampleFromOther(TWindow & other, TSiteSubset & subset, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);

	//getters
	GenotypeLikelihoods::TSite& operator[](uint32_t internalPos){ return _sites[internalPos]; };
	const std::string& chrName() const{ return _chrName; };
	uint32_t refId() const{ return _refId; };
	uint32_t posInRef(uint32_t internalPos) const{ return startPos + internalPos; };
	double depth();
	double fractionSitesNoData();
	double fractionDepthAtLeastTwo();
	uint32_t numSitesWithData();
	double fractionRefIsN();
	void dataSummary(TLog* Logfile);
	bool filter(const double maxFracMissing, const double maxRefN, TLog* Logfile);
	bool passedFilters() const{ return _passedFilters; };

	BAM::TAlignment* swapUsedForEmptyAlignment(BAM::TAlignment* usedAlignment);
	void initSites(long newLength);
	void clear();
	virtual void move(unsigned int Start, unsigned int End, int chrNumber, TLog* logfile);
	void addReferenceBaseToSites(BAM::TFastaBuffer & reference);
	void addReferenceBaseToSites(TSiteSubset & subset);
	void applyMask(BAM::TBedReaderWindows* mask, bool inverseMasking);
	void maskCpG(BAM::TFastaBuffer & reference);
	void estimateBaseFrequencies(GenotypeLikelihoods::TBaseData & baseFreq);
	void applyDepthFilter(const size_t minDepth, const size_t maxDepth);

	//loop over sites
	std::vector<GenotypeLikelihoods::TSite>::iterator begin(){ return _sites.begin(); };
	std::vector<GenotypeLikelihoods::TSite>::iterator end(){ return _sites.end(); };

	void writeCoordinates(TOutputFile & out);
};

//---------------------------------------------------------------
//TWindow
//---------------------------------------------------------------
class TWindow:public TWindow_base{
friend class TWindow_base;
private:
	//alignment stacks
	std::vector<BAM::TAlignment*> usedAlignments;
	std::vector<BAM::TAlignment*> emptyAlignments;

	void _cleanUpUsedAlignments(TLog* logfile);
	void _clearAllUsedAlignments();

	//functions to fill sites from alignments
	void _checkAlignmentForFillingSites(BAM::TAlignment* alignmentIt);
	void _setFirstPositionWithinWindow(BAM::TAlignment* alignmentIt, uint16_t & p);

	void _fillSites(BAM::TAlignment* alignmentIt, std::vector<GenotypeLikelihoods::TSite> & sites, const unsigned int & readUpToDepth);
	void _fillSitesSubset(BAM::TAlignment* alignmentIt, std::vector<GenotypeLikelihoods::TSite> & sites, std::set<TSiteSubsetSite> & thesePos, const unsigned int & readUpToDepth);

	void _fillSites(std::vector<GenotypeLikelihoods::TSite> & sites, const uint32_t & readUpToDepth);
	int _fillSitesDownsampling(std::vector<GenotypeLikelihoods::TSite> & sites, const uint32_t & readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);

	void _fillSitesSubset(std::vector<GenotypeLikelihoods::TSite> & sites, TSiteSubset & subset, const uint32_t & readUpToDepth);
	int _fillSitesSubsetDownsampling(std::vector<GenotypeLikelihoods::TSite> & sites, TSiteSubset & subset, const uint32_t & readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);

public:
	TWindow();
	~TWindow();

	void move(unsigned int Start, unsigned int End, int RefId, TLog* logfile);
	void review();
	void printStacks();

	void fillSites(const uint32_t & readUpToDepth);
	void fillSitesSubset(TSiteSubset & subset, const uint32_t & readUpToDepth);

	BAM::TAlignment* swapUsedForEmptyAlignment(BAM::TAlignment* usedAlignment);
};

}; //end namespace

#endif /* TWINDOW_H_ */
