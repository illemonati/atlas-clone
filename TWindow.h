/*
 * TWindow.h
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#ifndef TWINDOW_H_
#define TWINDOW_H_

#include "TLog.h"
#include "TParameters.h"
#include "TReadGroups.h"
#include "TThetaEstimator.h"
#include "TBedReader.h"
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


	std::string _chrName;

	TGenotypeMap genoMap;
	bool referenceBaseAdded;

	void setCoordinates(long Start, long End, int ChrNumber);

	//TODO: make as much as possible private
public:
	unsigned int startPos;
	unsigned int endPos; //end NOT included in window!
	unsigned int length;
	unsigned int chrNumber;

	std::vector<TSite> sites;
	bool sitesInitialized;
	unsigned int numReadsInWindow;
	double depth, fractionSitesNoData, fractionDepthAtLeastTwo;
	double fractionRefIsN;
	unsigned int numSitesWithData;
	bool passedFilters;
	TBaseFrequencies baseFreq;

	TWindow_base();
	TWindow_base(TWindow_base & other);
	void stealFromOther(TWindow_base & other);
	TWindow_base(TWindow & other, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);
	void downsampleFromOther(TWindow & other, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);
	void downsampleFromOther(TWindow & other, TSiteSubset & subset, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);

	virtual ~TWindow_base();

	//getters
	const std::string& chrName() const{ return _chrName; };
	uint32_t posInRef(uint32_t internalPos) const{ return startPos + internalPos; };
	TBaseFrequencies getBaseFreq(){return baseFreq;};


	BAM::TAlignment* swapUsedForEmptyAlignment(BAM::TAlignment* usedAlignment);
	void initSites(long newLength);
	void clear();
	virtual void move(unsigned int Start, unsigned int End, int chrNumber, TLog* logfile);
	void addReferenceBaseToSites(BAM::TFastaBuffer & reference);
	void addReferenceBaseToSites(TSiteSubset & subset);
	void applyMask(BAM::TBedReader* mask, bool inverseMasking);
	void maskCpG();
	void estimateBaseFrequencies();
	void calcDepth();
	void calcFracN();
	void calcDepthPerSite(long * siteDepth, size_t maxCov);
	void countDepthPerSite(TDistributionOfCounts & counts);
	void printDepthPerSite(gz::ogzstream & out, const std::string & chr);
	void countAlleles(TAllelicDepthCounts & counts);
	void writeNonConservedBed(std::ofstream & output);
	void applyDepthFilter(const size_t minDepth, const size_t maxDepth);
	void createDepthMask(size_t minDepth, size_t maxDepth, std::ofstream & outputMaskFile, const std::string & chr);

	//loop over sites
	std::vector<TSite>::iterator begin(){ return sites.begin(); };
	std::vector<TSite>::iterator end(){ return sites.end(); };

	//add sites to data structures
	void addSitesToThetaEstimator(TThetaEstimatorData* thetaDataContainer, TGenotypeLikelihoodCalculator & GL_calculator);
	void addSitesToThetaEstimator(TThetaEstimatorData* thetaDataContainer, TGenotypeLikelihoodCalculator & GL_calculator, BAM::TBedReader & region);
	void addToGLF(TGlfWriter & writer, TGenotypeLikelihoodCalculator & GL_calculator, bool printAll);
	void addToRecalibrationEM(TRecalibrationEMEstimator & recalObject, TQualityMap & qualMap);
	void addToRecalibrationEM(TRecalibrationEMEstimator & recalObject, TSiteSubset & subset, TQualityMap & qualMap);

	//other
	void generatePSMCInput(TThetaEstimator & estimator, TGenotypeLikelihoodCalculator & GL_calculator, int & blockSize, double & confidence, std::ofstream & out, int & nCharOnLine);
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
