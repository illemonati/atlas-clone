/*
 * TGenome_windows.h
 *
 *  Created on: May 31, 2020
 *      Author: phaentu
 */

#ifndef TGENOMEWINDOWS_H_
#define TGENOMEWINDOWS_H_

#include "TGenome.h"


//---------------------------------------------------------------
//TGenomeWindows
//---------------------------------------------------------------
class TGenomeWindows:public TGenome_recalibrated{
private:
	BAM::TBamFile bamFile;
 	TAlignmentParser alignmentParser;

	BAM::TFastaBuffer reference;




	void jumpToEnd();

	TGenotypePrior* initializeGenotypePrior(TParameters & params);
	void openSiteSubset(TBedReader* subset, std::string filename);
	bool estimateTheta(TThetaEstimator & thetaEstimator, TWindow_base & window);
	void fillVectorOfDownsamplingProbabilities(std::string prob, std::vector<double> & downSampleProbVector);
	void renameBAMSIfDouble(std::vector<double> & fracVector, std::vector<std::string> & names);

public:
	TGenomeWindows(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	~TGenomeWindows(){};

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

	//other
	void writeGLF(TParameters & params);
	void printPileup(TParameters & params);

	void estimatePMD(TParameters & params);
	void runPMDS(TParameters & params);
	void generatePSMCInput(TParameters & params);

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





#endif /* TGENOMEWINDOWS_H_ */
