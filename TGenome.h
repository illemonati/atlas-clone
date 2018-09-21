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
#include <typeinfo>
#include <map>
#include "TAlignmentParser.h"

//---------------------------------------------------------------
//TGenome
//---------------------------------------------------------------
class TGenome{
private:
 	TAlignmentParser alignmentParser;
/*
	TPMD* pmdObjects;
	bool hasPMD;
	TRecalibration* recalObject;
	TRecalibration* recalObject2;
	bool doRecalibration;
	bool doRecalibration2;
	bool recalObjectInitialized;
	bool recalObjectInitialized2;
*/
//	BamTools::Fasta reference;
//	bool fastaReference;
 	TRandomGenerator* randomGenerator;
 	bool randomGeneratorInitialized;

	TLog* logfile;
//	long oldPos;
	std::string outputName;
	int maxReadLength;

	void jumpToEnd();

	void openThetaOutputFile(std::ofstream & out, TThetaEstimator & estimator);
	void initializeRandomGenerator(TParameters & params);
	void openSiteSubset(TBedReader* subset, std::string filename);
	void indexBamFile(std::string & filename);

public:
	TGenome(TLog* Logfile, TParameters & params);
	~TGenome(){
		if(fastaReference) reference.Close();
//		if(recalObjectInitialized) delete recalObject;
//		if(pmdObjects) delete[] pmdObjects;
		if(randomGeneratorInitialized) delete randomGenerator;
	};

	//theta estimation
	bool initThetaEstimatorForCallers(TParameters & params, TThetaEstimator* & thetaEstimator);
	void estimateTheta(TParameters & params);
	void estimateThetaWindows(TThetaEstimator & thetaEstimator, std::ofstream & out);
	void estimateThetaGenomeWide(TThetaEstimator & thetaEstimator, std::ofstream & out, bool onlyReadData);
	void bootstrapTetaEstimation(int numBootstraps, TThetaEstimator & thetaEstimator);
	void calcLikelihoodSurfaces(TParameters & params);

	//callers
	bool openFastaReferenceForCaller(TParameters & params, BamTools::Fasta & reference);
	void callMLEGenotypes(TParameters & params);
	void callBayesianGenotypes(TParameters & params);
	void callAllelePresence(TParameters & params);
	void randomBaseCaller(TParameters & params);
	void majorityBaseCaller(TParameters & params);

	//other
	void writeGLF(TParameters & params);
	void combineBeagleFiles(TParameters & params);
	void printPileup(TParameters & params);

	//recalibration
	void estimateErrorCalibration(TParameters & params);
	void estimateErrorCalibrationEM(TParameters & params);
	void calculateLikelihoodErrorCalibrationEM(TParameters & params);
	void BQSR(TParameters & params);
	void printQualityDistribution(TParameters & params);
	void printQualityTransformation(TParameters & params);
	void reportProgressParsingBamFile(const long & counter, const struct timeval & start);
	void recalibrateBamFile(TParameters & params);
	void binQualityScores(TParameters & params);
	void assessSoftClipping(TParameters & params);
	void assessOverlap(TParameters & params);
	void splitSingleEndReadGroups(TParameters & params);
	void mergeReadGroups(TParameters & params);
	void estimatePMD(TParameters & params);
	float calculatePMDS(int readGroup, char & ref, char & read, double & pmdCT, double & pmdGA, double & errorRate, double & pi, float & probPMD, float & probNoPMD);
	void runPMDS(TParameters & params);
	void mergePairedEndReads(TParameters & params);
	void generatePSMCInput(TParameters & params);
	void downSampleBamFile(TParameters & params);
	void downSampleReads(TParameters & params);
	void diagnoseBamFile(TParameters & params);
	void allelicDepth(TParameters & params);
	void estimateApproximateDepthPerWindow(TParameters & params);
	void estimateDepthPerSite(TParameters & params);
	void writeDepthPerSite(TParameters & params);
	void createDepthMask(TParameters & params);
	void simulateGWASData(TParameters & params);
};




#endif /* LOCI_H_ */
