/*
 * loci.h
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#ifndef LOCI_H_
#define LOCI_H_

#include "TWindow.h"
#include "TRecalibration.h"
#include "gzstream.h"
#include "bamtools/api/BamWriter.h"
#include "TLog.h"
#include "TBed.h"
#include "TAlignmentParser.h"
#include "TQualityMap.h"
#include "TGenotypePrior.h"
#include <typeinfo>
#include <map>
#include <algorithm>

//---------------------------------------------------------------
//TGenome
//---------------------------------------------------------------
class TGenome{
private:
 	TAlignmentParser alignmentParser;
	BamTools::Fasta reference;
 	TRandomGenerator* randomGenerator;
 	bool randomGeneratorInitialized;

	TLog* logfile;
	std::string outputName;
	int maxReadLength;

	void jumpToEnd();

	void initializeRandomGenerator(TParameters & params);
	TGenotypePrior* initializeGenotypePrior(TParameters & params);
	void openSiteSubset(TBedReader* subset, std::string filename);
	void indexBamFile(std::string & filename);
	void mergeAlignedBasesBamReads(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality);
	void dealWithLastReadsInStorage(std::vector< std::pair<TAlignment*, bool> > & alignmentStorage, const bool & filterOrphanedReads);

public:
	TGenome(TLog* Logfile, TParameters & params);
	~TGenome(){
		if(randomGeneratorInitialized) delete randomGenerator;
	};

	//theta estimation
	bool initThetaEstimatorForCallers(TParameters & params, TThetaEstimator* & thetaEstimator);
	void estimateTheta(TParameters & params);
	void estimateThetaWindows(TThetaEstimator & thetaEstimator, TThetaOutputFile & out);
	void estimateThetaGenomeWide(TThetaEstimator & thetaEstimator, TThetaOutputFile & out, bool onlyReadData, int numBootstraps);
	void bootstrapTetaEstimation(int numBootstraps, TThetaEstimator & thetaEstimator);
	void calcThetaLikelihoodSurfaces(TParameters & params);
	void estimateThetaRatio(TParameters & params);

	//callers
	void callGenotypes(TParameters & params);

//	bool openFastaReferenceForCaller(TParameters & params, BamTools::Fasta & reference);
//	void writeVcfHeader(gz::ogzstream* output, bool limitToSitesWithKnownAlleles, bool onlyPhredGP);
//	void callMLEGenotypes(TParameters & params);
//	void callBayesianGenotypes(TParameters & params);
//	void callAllelePresence(TParameters & params);
//	void randomBaseCaller(TParameters & params);
//	void majorityBaseCaller(TParameters & params);


	//recalibration
	void estimateErrorCalibration(TParameters & params);
	void estimateErrorCalibrationEM(TParameters & params);
	void calculateLikelihoodErrorCalibrationEM(TParameters & params);
	void BQSR(TParameters & params);
	void printQualityDistribution(TParameters & params);
	void printQualityTransformation(TParameters & params);

	//other
	void writeGLF(TParameters & params);
	void printPileup(TParameters & params);
	void reportProgressParsingBamFile(const long & counter, const struct timeval & start);
	void reportProgressParsingBamFileNoCheck(const long & counter, const struct timeval & start);
	void recalibrateBamFile(TParameters & params);
	void binQualityScores(TParameters & params);
	void assessSoftClipping(TParameters & params);
	void removeSoftClippedBasesFromReads(TParameters & params);
	void assessOverlap(TParameters & params);
	void splitSingleEndReadGroups(TParameters & params);
	void mergeReadGroups(TParameters & params);
	void estimatePMD(TParameters & params);
	float calculatePMDS(int readGroup, char & ref, char & read, double & pmdCT, double & pmdGA, double & errorRate, double & pi, float & probPMD, float & probNoPMD);
	void runPMDS(TParameters & params);
	void mergePairedEndReads(TParameters & params);
	void updateOrphanedReadsAtBeginningOfStorage(std::vector< std::pair<TAlignment*, bool> > & alignmentStorage, TAlignment & alignment, int & acceptedDistanceBetweenMates, const bool & filterOrphanedReads);
	void writeAllReadsThatAreReady(BamTools::BamWriter & bamWriter, std::vector< std::pair<TAlignment*, bool> > & alignmentStorage);
	void findPairedReadGroupsToMergeReads(TParameters & params, std::vector<bool> & pairedReadGroups, bool & allReadGroupsPaired);
	bool ignoreReadAfterChrSwitch(std::vector< std::pair<TAlignment*, bool> > & alignmentStorage, TAlignment & alignment, const bool & filterOrphanedReads);
	bool findAndMergeMates(std::vector< std::pair<TAlignment*, bool> > & alignmentStorage, TAlignment & alignment, const bool & adaptQuality);
	void decideWhatToDoWithLastReadsInStorage(std::vector< std::pair<TAlignment*, bool> > & alignmentStorage, const bool & filterOrphanedReads);
	void mergePairedEndReadsNoOrder(TParameters & params);
	void generatePSMCInput(TParameters & params);
	void downSampleBamFile(TParameters & params);
	void downSampleReads(TParameters & params);
	void diagnoseBamFile(TParameters & params);
	void allelicDepth(TParameters & params);
	void estimateApproximateDepthPerWindow(TParameters & params);
	void estimateDepthPerSite(TParameters & params);
	void writeDepthPerSite(TParameters & params);
	void estimateDuplicationCounts(TParameters & params);
	void createDepthMask(TParameters & params);
	void simulateGWASData(TParameters & params);
	void printMateInformationPerSite(TParameters & params);
};




#endif /* LOCI_H_ */
