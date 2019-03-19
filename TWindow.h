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
#include "TRecalibration.h"
#include "TRecalibrationBQSR.h"
#include "TThetaEstimator.h"
#include "TBedReader.h"
#include "TSiteSubset.h"
#include "TPostMortemDamage.h"
#include "TGLF.h"
#include "TAlignment.h"
#include "TQualityMap.h"
#include "TCaller.h"
#include "TDistributionOfCounts.h"
#include "TRecalibrationEMEstimator.h"


//---------------------------------------------------------------
//TWindow
//---------------------------------------------------------------
class TWindow{
private:
	TGenotypeMap genoMap;
	bool referenceBaseAdded;

	//alignment stacks
	std::vector<TAlignment*> usedAlignments;
	std::vector<TAlignment*> emptyAlignments;

	void fillPGenotype(double* pGenotype);
	//std::vector<TAlignment*>::iterator lastAlignmentwithEndInWindow;
	//std::vector<TAlignment*>::iterator firstAlignmentwithPosOutsideWindow;

public:
	long start;
	long end; //end NOT included in window!
	long length;
	int chrNumber;
	std::string chrName;
	TSite* sites;
	bool sitesInitialized;
	int numReadsInWindow;
	double depth, fractionSitesNoData, fractionsitesDepthAtLeastTwo;
	double fractionRefIsN;
	long numSitesWithData;
	bool passedFilters;
	TBaseFrequencies baseFreq;

	TWindow();
	~TWindow();

	//getters
	TBaseFrequencies getBaseFreq(){return baseFreq;};

	TAlignment* swapUsedForEmptyAlignment(TAlignment* usedAlignment, const unsigned int & maxReadLength);
	void initSites(long newLength);
	void clear();
	void move(long Start, long End, int chrNumber);
	void review();
	void cleanUpUsedAlignments();
	void printStacks();
	void fillSitesSubset(TSiteSubset* subset, const int & readUpToDepth);
	void fillSites(const int & readUpToDepth);
	void addReferenceBaseToSites(BamTools::Fasta & reference, int & refId);
	void addReferenceBaseToSites(TSiteSubset* subset);
	void applyMask(TBedReader* mask, bool inverseMasking);
	void maskCpG(BamTools::Fasta & reference, int & refId);
	void estimateBaseFrequencies();
	void calculateEmissionProbabilities();
	void callMLEGenotype(TRecalibration* recalObject, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF, bool gVCF, bool noAltIfHomoRef);
	void printPileup(TRecalibration* recalObject, gz::ogzstream & out, bool printOnlySitesWithData);
	void printPileupToScreen(TRecalibration* recalObject);
	void calcDepth();
	void calcFracN();
	void calcDepthPerSite(long * siteDepth, size_t maxCov);
	void countDepthPerSite(TDistributionOfCounts & counts);
	void printDepthPerSite(gz::ogzstream & out, std::string & chr);
	void countAlleles(long**** siteImbalance, const unsigned int & maxDepth);
	void applyDepthFilter(int minDepth, size_t maxDepth);
	void createDepthMask(size_t minDepth, size_t maxDepth, std::ofstream & outputMaskFile, std::string & chr);

	//add sites to data structures
	void addSitesToBQSR(TRecalibrationBQSREstimator & bqsr, TLog* logfile, TQualityMap & qualMap);
	void addSitesToBQSR(TRecalibrationBQSREstimator & bqsr, TSiteSubset* subset, TLog* logfile, TQualityMap & qualMap);
	void addSitesToQualityTransformTable(TRecalibration* recalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile, TQualityMap & qualMap);
	void addSitesToQualityTransformTable(TRecalibration* recalObject, TRecalibration* otherRecalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile, TQualityMap & qualMap);
	void addSitesToPMDTable(TPMDTables & pmdTables, TLog* logfile);
	void addSitesToThetaEstimator(TThetaEstimatorData* thetaDataContainer);
	void addSitesToThetaEstimator(TThetaEstimatorData* thetaDataContainer, TBedReader & region);
	void addToGLF(TGlfWriter & writer, bool printAll);
	void addToRecalibrationEM(TRecalibrationEMEstimator & recalObject, TQualityMap & qualMap);
	void addToRecalibrationEM(TRecalibrationEMEstimator & recalObject, TSiteSubset* subset, TQualityMap & qualMap);

	//callers
	void call(TCaller & caller, TRecalibration & recalObject, BamTools::Fasta & reference);
	void callKnwonAlleles(TCaller & caller, TRecalibration & recalObject, BamTools::Fasta & reference, TSiteSubset & subset);

	//other
	void generatePSMCInput(TThetaEstimator & estimator, int & blockSize, double & confidence, std::ofstream & out, int & nCharOnLine);
	double calcLogLikelihood();
};

#endif /* TWINDOW_H_ */
