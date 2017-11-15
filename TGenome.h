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

//---------------------------------------------------------------
//TGenome
//---------------------------------------------------------------
class TGenome{
private:
	TPMD* pmdObjects;
	bool hasPMD;
	TRecalibration* recalObject;
	TRecalibration* recalObject2;
	bool doRecalibration;
	bool doRecalibration2;
	bool recalObjectInitialized;
	bool recalObjectInitialized2;

	BamTools::BamReader bamReader;
	BamTools::BamRegion bamRegion;
 	BamTools::SamHeader bamHeader;
 	BamTools::BamAlignment bamAlignment;
 	TAlignmentParser alignmentParser;
	bool oldAlignementMustBeConsidered;
	BamTools::Fasta reference;
	bool fastaReference;
 	BamTools::SamSequenceIterator chrIterator;
 	TReadGroups readGroups;
 	TRandomGenerator* randomGenerator;
 	bool randomGeneratorInitialized;
 	int chrNumber;
 	long chrLength;
 	long curStart;
 	long curEnd;
	std::string filename;
	TLog* logfile;
	bool windowsPredefined;
	TBed* predefinedWindows;
	int windowSize;
	int numWindowsOnChr;
	int windowNumber;
	int maxReadLength;
	double maxMissing;
	double maxRefN;
	long oldPos;
	std::string outputName;
	TBedReader* mask;
	bool doMasking, considerRegions;
	bool doCpGMasking;
	bool applyCoverageFilter, applyQualityFilter;
	size_t minCoverage, maxCoverage;
	int minQuality, maxQuality;
	int minOutQuality, maxOutQuality;
	long limitWindows;
	int limitChr;
	bool* useChromosome;
	bool limitReadGroups;
	std::vector<std::string> readGroupsInUse;

	void jumpToEnd();
	void restartChromosome(TWindowPair & windowPair);
	bool iterateChromosome(TWindowPair & windowPair);
	void moveChromosome(TWindowPair & windowPair);
	bool iterateWindow(TWindowPair & windowPair);
	bool addAlignementToWindows(TAlignmentParser & alignment, TWindowPair & windowPair);
	bool readData(TWindowPair & windowPair);
	void initializePostMortemDamage(TParameters & params);
	void initializeRecalibration(TParameters & params);
	void openThetaOutputFile(std::ofstream & out, TThetaEstimator & estimator);
	void initializeRandomGenerator(TParameters & params);
	void openSiteSubset(TBedReader* subset, std::string filename);

public:
	TGenome(TLog* Logfile, TParameters & params);
	~TGenome(){
		if(doMasking){
			std::cout << "----------------------__> DELETE!!!!" << std::endl;
			delete mask;
		}
		if(fastaReference) reference.Close();
		if(recalObjectInitialized) delete recalObject;
		if(pmdObjects) delete[] pmdObjects;
		if(randomGeneratorInitialized) delete randomGenerator;
		if(useChromosome) delete[] useChromosome;
		if(windowsPredefined) delete predefinedWindows;
	};

	//theta estimation
	bool initThetaEstimatorForCallers(TParameters & params, TThetaEstimator* & thetaEstimator);
	void estimateTheta(TParameters & params);
	void estimateThetaWindows(TThetaEstimator & thetaEstimator, std::ofstream & out);
	void estimateThetaGenomeWide(TThetaEstimator & thetaEstimator, std::ofstream & out);
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
	void estimateErrorCalibration(TParameters & params);
	void estimateErrorCalibrationEM(TParameters & params);
	//void fillSequence(std::vector<double> & vec, std::string & str);
	void calculateLikelihoodErrorCalibrationEM(TParameters & params);
	void BQSR(TParameters & params);
	void printQualityTransformation(TParameters & params);
	void createBase(TBase** basePointer, char & base, char & quality, int & posInRead, int & revPosInRead, double & pmdCT, double & pmdGA, BaseContext & context, int & readGroupId);
	char returnBaseQualityAsChar(char & base, char & quality, int & posInRead, int & revPosInRead, double & pmdCT, double & pmdGA, BaseContext & context, int & readGroupId);
	double returnBaseQualityWithPMDAsCharFwdMapping(char & base, char & refBase, char & quality, int & posInRead, int & revPosInRead, double & pmdCT, double & pmdGA, BaseContext & context, int & readGroupId);
	double returnBaseQuality(char & base, char & quality, int & posInRead, int & revPosInRead, double & pmdCT, double & pmdGA, BaseContext & context, int & readGroupId);
	double returnBaseQualityWithPMDAsCharRevMapping(char & base, char & refBase, char & quality, int & posInRead, int & revPosInRead, double & pmdCT, double & pmdGA, BaseContext & context, int & readGroupId);
	bool recalibrateAlignment(BamTools::BamAlignment & alignment, std::string & qual, TGenotypeMap & genoMap, bool withPMD, int & begin, std::string & ref, std::map <std::string, int> & mateTooLong);
	void recalibrateBamFile(TParameters & params);
	void binQualityScores(TParameters & params);
	void assessSoftClipping(TParameters & params);
	void splitSingleEndReadGroups(TParameters & params);
	void mergeReadGroups(TParameters & params);
	void addReadToPMD(TWindowDiploid* window, TGenotypeMap & genoMap, std::string & ref, TPMDTables & pmdTables);
	void estimatePMD(TParameters & params);
	float calculatePMDS(int readGroup, char & ref, char & read, double & pmdCT, double & pmdGA, double & errorRate, double & pi, float & probPMD, float & probNoPMD);
	void runPMDS(TParameters & params);
	void mergePairedEndReads(TParameters & params);
	void generatePSMCInput(TParameters & params);
	void downSampleBamFile(TParameters & params);
	void downSampleReads(TParameters & params);
	void diagnoseBamFile(TParameters & params);
	void estimateApproximateCoveragePerWindow(TParameters & params);
	void estimateCoveragePerSite(TParameters & params);
	void createDepthMask(TParameters & params);
	void simulateGWASData(TParameters & params);
};




#endif /* LOCI_H_ */
