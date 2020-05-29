/*
 * loci.h
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#ifndef LOCI_H_
#define LOCI_H_

#include <TGenotypeData.h>

#include "TWindow.h"
#include "gzstream.h"
#include "bamtools/api/BamWriter.h"
#include "TLog.h"
#include "TBed.h"
#include "TAlignmentParser.h"
#include "TQualityMap.h"
#include "TReadList.h"
#include <typeinfo>
#include <map>
#include <algorithm>
#include "counters.h"

#include "TRecalibration.h"
#include "TAlignmentMerger.h"
#include "TSoftClipping.h"
#include "TAllelicDepthCounts.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TGenotypePrior.h"

//---------------------------------------------------------------
//TGenome
//---------------------------------------------------------------
class TGenome{
private:
	BAM::TBamFile bamFile;
 	TAlignmentParser alignmentParser;
 	TGenotypeLikelihoodCalculator genotypeLikelihoodCalculator;
	BAM::TFastaBuffer reference;
 	TRandomGenerator* randomGenerator;

	TLog* logfile;
	std::string outputName;

	void jumpToEnd();

	TGenotypePrior* initializeGenotypePrior(TParameters & params);
	void openSiteSubset(TBedReader* subset, std::string filename);
	void indexBamFile(std::string & filename);
	bool estimateTheta(TThetaEstimator & thetaEstimator, TWindow_base & window);
	void mergeAlignedBasesBamReads(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality);
	void findPairedReadGroupsToMergeReads(TParameters & params, std::vector<bool> & pairedReadGroups);
	void parseSplitMergeReadGroupSettings(TParameters & params, std::map<int, TReadGroupMaxLength> & RGSettings);
	void setMergerSettings(TParameters & params, TAlignmentMerger & merger);
	void fillVectorOfDownsamplingProbabilities(std::string prob, std::vector<double> & downSampleProbVector);
	void renameBAMSIfDouble(std::vector<double> & fracVector, std::vector<std::string> & names);

public:
	TGenome(TLog* Logfile, TParameters & params, TRandomGenerator* RandomGenerator);
	~TGenome(){};

	//theta estimation
	bool initThetaEstimatorForCallers(TParameters & params, TThetaEstimator* & thetaEstimator);
	void estimateTheta(TParameters & params);
	void estimateThetaWindows(TThetaEstimator & thetaEstimator, TThetaOutputFile & out, bool printAll);
	void estimateThetaGenomeWide(TThetaEstimator & thetaEstimator, TThetaOutputFile & out, bool onlyReadData, int numBootstraps);
	void bootstrapTetaEstimation(int numBootstraps, TThetaEstimator & thetaEstimator);
	void calcThetaLikelihoodSurfaces(TParameters & params);
	void estimateThetaRatio(TParameters & params);
	void performDownsamplingThetaQC(TParameters & params);

	//callers
	void callGenotypes(TParameters & params);

	//recalibration
	void estimateErrorCalibrationEM(TParameters & params);
	void calculateLikelihoodErrorCalibrationEM(TParameters & params);
	void BQSR(TParameters & params);
	void printQualityDistribution(TParameters & params);
	void printQualityTransformation(TParameters & params);

	//other
	void writeGLF(TParameters & params);
	void printPileup(TParameters & params);
	void recalibrateBamFile(TParameters & params);
	void binQualityScores(TParameters & params);
	void assessSoftClipping(TParameters & params);
	void removeSoftClippedBasesFromReads(TParameters & params);
	void assessOverlap(TParameters & params);
	void splitSingleEndReadGroups(TParameters & params);
	void mergeReadGroups(TParameters & params);
	void estimatePMD(TParameters & params);
	void runPMDS(TParameters & params);
	void filterBAM(TParameters & params);
	void splitMerge(TParameters & params);
	void mergePairedEndReads(TParameters & params);
	void generatePSMCInput(TParameters & params);
	void downSampleBamFile(TParameters & params);
	void separateReads(TParameters & params);
	void downSampleReads(TParameters & params);
	void diagnoseBamFile(TParameters & params);
	void allelicDepth(TParameters & params);
	void writeNonConservedBed(TParameters & params);
	void estimateApproximateDepthPerWindow(TParameters & params);
	void estimateDepthPerSite(TParameters & params);
	void writeDepthPerSite(TParameters & params);
	void estimateDuplicationCounts(TParameters & params);
	void createDepthMask(TParameters & params);
	void printMateInformationPerSite(TParameters & params);
	void contextStats(TParameters & params);

	//debug stuff
	void testGenotypeLikelihoods(TParameters & params);
};




#endif /* LOCI_H_ */
