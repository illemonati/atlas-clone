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
#include "TCaller.h"
#include "TDistributionOfCounts.h"
#include "TRandomGenerator.h"
#include "GenotypeLikelihoods/TRecalibration.h"
#include "GenotypeLikelihoods/TRecalibrationBQSR.h"
#include "GenotypeLikelihoods/TRecalibrationEMEstimator.h"

//forward declaration to enable copy constructor
class TWindow;

//---------------------------------------------------------------
//TWindow_base
//---------------------------------------------------------------
class TWindow_base{
protected:
	TGenotypeMap genoMap;
	bool referenceBaseAdded;

	void setCoordinates(long Start, long End, int ChrNumber);

public:
	unsigned int start;
	unsigned int end; //end NOT included in window!
	unsigned int length;
	unsigned int chrNumber;
	std::string chrName;
	TSite* sites;
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
	void downsampleFromOther(TWindow & other, TSiteSubset* subset, const int readUpToDepth, const double downsamplingProb, TRandomGenerator* randomGenerator);

	virtual ~TWindow_base();

	//getters
	TBaseFrequencies getBaseFreq(){return baseFreq;};

	TAlignment* swapUsedForEmptyAlignment(TAlignment* usedAlignment, const unsigned int & maxReadLength);
	void initSites(long newLength);
	void clear();
	virtual void move(unsigned int Start, unsigned int End, int chrNumber);
	void addReferenceBaseToSites(BamTools::Fasta & reference);
	void addReferenceBaseToSites(TSiteSubset* subset);
	void applyMask(TBedReader* mask, bool inverseMasking);
	void maskCpG();
	void estimateBaseFrequencies();
	void calculateEmissionProbabilities();
	void callMLEGenotype(TRecalibration* recalObject, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF, bool gVCF, bool noAltIfHomoRef);
	void printPileup(TRecalibration* recalObject, gz::ogzstream & out, bool printOnlySitesWithData);
	void printPileupToScreen(TRecalibration* recalObject);
	void calcDepth();
	void calcFracN();
	void calcDepthPerSite(long * siteDepth, size_t maxCov);
	void countDepthPerSite(TDistributionOfCounts & counts);
	void printDepthPerSite(gz::ogzstream & out, const std::string & chr);
	void printMateInformationPerSite(TOutputFileZipped & out, const std::string & chr);
	void countAlleles(TAllelicDepthCounts & counts);
	void writeNonConservedBed(std::ofstream & output);
	void contextStats(int** contextCounts, TQualityMap & qualMap);
	void applyDepthFilter(const size_t minDepth, const size_t maxDepth);
	void createDepthMask(size_t minDepth, size_t maxDepth, std::ofstream & outputMaskFile, const std::string & chr);

	//add sites to data structures
	void addSitesToBQSR(TRecalibrationBQSREstimator & bqsr, TLog* logfile, TQualityMap & qualMap);
	void addSitesToBQSR(TRecalibrationBQSREstimator & bqsr, TSiteSubset* subset, TLog* logfile, TQualityMap & qualMap);
	void addSitesToQualityTransformTable(TRecalibration* recalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile, TQualityMap & qualMap);
	void addSitesToQualityTransformTable(TRecalibration* recalObject, TRecalibration* otherRecalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile, TQualityMap & qualMap);
	void addSitesToPMDTable(GenotypeLikelihoods::TPMDTables & pmdTables, TLog* logfile);
	void addSitesToThetaEstimator(TThetaEstimatorData* thetaDataContainer);
	void addSitesToThetaEstimator(TThetaEstimatorData* thetaDataContainer, TBedReader & region);
	void addToGLF(TGlfWriter & writer, const int ploidy, bool printAll);
	void addToRecalibrationEM(GenotypeLikelihoods::TRecalibrationEMEstimator & recalObject, TQualityMap & qualMap);
	void addToRecalibrationEM(GenotypeLikelihoods::TRecalibrationEMEstimator & recalObject, TSiteSubset* subset, TQualityMap & qualMap);

	//callers
	void call(TCaller & caller, TRecalibration & recalObject, BamTools::Fasta & reference);
	void callKnwonAlleles(TCaller & caller, TRecalibration & recalObject, TSiteSubset & subset);

	//other
	void generatePSMCInput(TThetaEstimator & estimator, int & blockSize, double & confidence, std::ofstream & out, int & nCharOnLine);
};

//---------------------------------------------------------------
//TWindow
//---------------------------------------------------------------
class TWindow:public TWindow_base{
private:
	//alignment stacks
	std::vector<TAlignment*> usedAlignments;
	std::vector<TAlignment*> emptyAlignments;

	void cleanUpUsedAlignments();
	void clearAllUsedAlignments();

	//functions to fill sites from alignments
	void checkAlignmentForFillingSites(TAlignment* alignmentIt);
	void setFirstPositionWithinWindow(TAlignment* alignmentIt, unsigned int & firstPos, unsigned int & p);
	void fillSites(TAlignment* alignmentIt, TSite* theseSites, const unsigned int & readUpToDepth);
	void fillSitesSubset(TAlignment* alignmentIt, TSite* theseSites, std::map<unsigned int,std::pair<char,char> > & thesePos, const unsigned int & readUpToDepth);

public:
	TWindow();
	~TWindow();

	void move(unsigned int Start, unsigned int End, int chrNumber);
	void review();
	void printStacks();

	void fillSites(const unsigned int & readUpToDepth);
	int fillSites(TSite* theseSites, const unsigned int & readUpToDepth);
	int fillSitesDownsampling(TSite* theseSites, const unsigned int & readUpToDepth, double downsamplingProb, TRandomGenerator* randomGenerator);

	void fillSitesSubset(TSiteSubset* subset, const unsigned int & readUpToDepth);
	int fillSitesSubset(TSite* theseSites, TSiteSubset* subset, const unsigned int & readUpToDepth);
	int fillSitesSubsetDownsampling(TSite* theseSites, TSiteSubset* subset, const unsigned int & readUpToDepth, double downsamplingProb, TRandomGenerator* randomGenerator);

	TAlignment* swapUsedForEmptyAlignment(TAlignment* usedAlignment, const unsigned int & maxReadLength);
};



#endif /* TWINDOW_H_ */
