/*
 * TWindow.h
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#ifndef TWINDOW_H_
#define TWINDOW_H_

#include <TBedReaderWindows.h>
#include "TLog.h"
#include "TParameters.h"
#include "TReadGroups.h"
#include "TThetaEstimator.h"
#include "TSiteSubset.h"
#include "TPostMortemDamage.h"
#include "GLF/TGLF.h"
#include "TAlignment.h"
#include "TQualityMap.h"
#include "TDistributionOfCounts.h"
#include "TRandomGenerator.h"
#include "GenotypeLikelihoods/TRecalibration.h"
#include "GenotypeLikelihoods/TRecalibrationEMEstimator.h"
#include "GenotypeLikelihoods/TGenotypeLikelihoodCalculator.h"

#include <vector>

//forward declaration to enable copy constructor
class TWindow;

//---------------------------------------------------------------
//TWindow_base
//---------------------------------------------------------------
class TWindow_base{
protected:
	std::vector<GenotypeLikelihoods::TSite> _sites;

	std::string _chrName;

	TGenotypeMap genoMap;
	bool referenceBaseAdded;

	//depth
	bool _depthCalculated;
	double _depth, _fractionSitesNoData, _fractionDepthAtLeastTwo, _fractionRefIsN;
	uint32_t _numSitesWithData;
	void _calcDepth();

	bool _passedFilters;

	void _setCoordinates(long Start, long End, int ChrNumber);

	//TODO: make as much as possible private
public:
	unsigned int startPos;
	unsigned int endPos; //end NOT included in window!
	unsigned int length;
	unsigned int chrNumber;


	bool sitesInitialized;
	unsigned int numReadsInWindow;

	TWindow_base();
	TWindow_base(TWindow_base & other);
	void stealFromOther(TWindow_base & other);
	TWindow_base(TWindow & other, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);
	void downsampleFromOther(TWindow & other, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);
	void downsampleFromOther(TWindow & other, TSiteSubset & subset, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);

	virtual ~TWindow_base();

	//getters
	GenotypeLikelihoods::TSite& operator[](uint32_t internalPos){ return _sites[internalPos]; };
	const std::string& chrName() const{ return _chrName; };
	uint32_t posInRef(uint32_t internalPos) const{ return startPos + internalPos; };
	double depth();
	double fractionSitesNoData();
	double fractionDepthAtLeastTwo();
	uint32_t numSitesWithData();
	double fractionRefIsN();
	void dataSummary(TLog* Logfile);
	bool filter(const double maxFracMissing, const double maxRefN, TLog* Logfile);
	bool passedFilters() const{ return _passedFilters; };

	void writeCoordinates(TOutputFile & out);

	BAM::TAlignment* swapUsedForEmptyAlignment(BAM::TAlignment* usedAlignment);
	void initSites(long newLength);
	void clear();
	virtual void move(unsigned int Start, unsigned int End, int chrNumber, TLog* logfile);
	void addReferenceBaseToSites(BAM::TFastaBuffer & reference);
	void addReferenceBaseToSites(TSiteSubset & subset);
	void applyMask(BAM::TBedReaderWindows* mask, bool inverseMasking);
	void maskCpG();
	void estimateBaseFrequencies(GenotypeLikelihoods::TBaseData & baseFreq);

	void countAlleles(TAllelicDepthCounts & counts);
	void writeNonConservedBed(std::ofstream & output);
	void applyDepthFilter(const size_t minDepth, const size_t maxDepth);

	//loop over sites
	std::vector<GenotypeLikelihoods::TSite>::iterator begin(){ return _sites.begin(); };
	std::vector<GenotypeLikelihoods::TSite>::iterator end(){ return _sites.end(); };

	//add sites to data structures
	void addToGLF(TGlfWriter & writer, TGenotypeLikelihoodCalculator & GL_calculator, bool printAll);
};

//---------------------------------------------------------------
//TWindow
//---------------------------------------------------------------
class TWindow:public TWindow_base{
private:
	//alignment stacks
	std::vector<BAM::TAlignment*> usedAlignments;
	std::vector<BAM::TAlignment*> emptyAlignments;

	void _cleanUpUsedAlignments(TLog* logfile);
	void _clearAllUsedAlignments();

	//functions to fill sites from alignments
	void _checkAlignmentForFillingSites(BAM::TAlignment* alignmentIt);
	void _setFirstPositionWithinWindow(BAM::TAlignment* alignmentIt, uint16_t & p);
	void _fillSites(BAM::TAlignment* alignmentIt, std::vector<TSite> & sites, const unsigned int & readUpToDepth);
	void _fillSitesSubset(BAM::TAlignment* alignmentIt, std::vector<TSite> & sites, std::set<TSiteSubsetSite> & thesePos, const unsigned int & readUpToDepth);

public:
	TWindow();
	~TWindow();

	void move(unsigned int Start, unsigned int End, int chrNumber, TLog* logfile);
	void review();
	void printStacks();

	void fillSites(const uint32_t & readUpToDepth);
	void fillSites(std::vector<TSite> & sites, const uint32_t & readUpToDepth);
	int fillSitesDownsampling(std::vector<TSite> & sites, const uint32_t & readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);

	void fillSitesSubset(TSiteSubset & subset, const uint32_t & readUpToDepth);
	void fillSitesSubset(std::vector<TSite> & sites, TSiteSubset & subset, const uint32_t & readUpToDepth);
	int fillSitesSubsetDownsampling(std::vector<TSite> & sites, TSiteSubset & subset, const uint32_t & readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);

	BAM::TAlignment* swapUsedForEmptyAlignment(BAM::TAlignment* usedAlignment);
};



#endif /* TWINDOW_H_ */
