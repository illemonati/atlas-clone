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

//---------------------------------------------------------------
//TGenome
//---------------------------------------------------------------
class TGenome{
private:
	void jumpToEnd();
	void restartChromosome(TWindowPair & windowPair);
	bool iterateChromosome(TWindowPair & windowPair);
	void moveChromosome(TWindowPair & windowPair);
	bool iterateWindow(TWindowPair & windowPair);
	bool addAlignementToWindows(BamTools::BamAlignment & alignement, TWindowPair & windowPair);
	bool readData(TWindowPair & windowPair);
	void initializePostMortemDamage(TParameters & params);
	void initializeRecalibration(TParameters & params);
	void openThetaOutputFile(std::ofstream & out);
	void initializeRandomGenerator(TParameters & params);

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
 	int windowSize;
 	int numWindowsOnChr;
 	int windowNumber;
 	double maxMissing;
 	long oldPos;
 	std::string outputName;
 	TBedReader* mask;
 	bool doMasking;
 	bool doCpGMasking;
 	bool applyCoverageFilter, applyQualityFilter;
 	int minCoverage, maxCoverage;
 	int minQuality, maxQuality;
 	long limitWindows;
 	int limitChr;
 	bool* useChromosome;

public:
	TGenome(TLog* Logfile, TParameters & params);
	~TGenome(){
		if(doMasking) delete mask;
		if(fastaReference) reference.Close();
		if(recalObjectInitialized) delete recalObject;
		if(pmdObjects) delete[] pmdObjects;
		if(randomGeneratorInitialized) delete randomGenerator;
		if(useChromosome) delete[] useChromosome;
	};
	void estimateTheta(TParameters & params);
	void calcLikelihoodSurfaces(TParameters & params);
	bool openFastaReferenceForCaller(TParameters & params, BamTools::Fasta & reference);
	void callMLEGenotypes(TParameters & params);
	void callBayesianGenotypes(TParameters & params);
	void callAllelePresence(TParameters & params);
	void printPileup(TParameters & params);
	void estimateErrorCalibration(TParameters & params);
	void estimateErrorCalibrationEM(TParameters & params);
	//void fillSequence(std::vector<double> & vec, std::string & str);
	void calculateLikelihoodSurfaceErrorCalibrationEM(TParameters & params);
	void BQSR(TParameters & params);
	void printQualityTransformation(TParameters & params);
	void createBase(TBase** basePointer, char & base, char & quality, int & posInRead, int & revPosInRead, double & pmdCT, double & pmdGA, BaseContext & context, int & readGroupId);
	char returnBaseQualityAsChar(char & base, char & quality, int & posInRead, int & revPosInRead, double & pmdCT, double & pmdGA, BaseContext & context, int & readGroupId);
	double returnBaseQuality(char & base, char & quality, int & posInRead, int & revPosInRead, double & pmdCT, double & pmdGA, BaseContext & context, int & readGroupId);
	bool recalibrateAlignment(BamTools::BamAlignment & alignment, std::string & qual, TGenotypeMap & genoMap, std::map <std::string, int> & mateTooLong);
	void recalibrateBamFile(TParameters & params);
	void splitSingleEndReadGroups(TParameters & params);
	void mergeReadGroups(TParameters & params);
	void addReadToPMD(TWindowDiploid* window, TGenotypeMap & genoMap, std::string & ref, TPMDTables & pmdTables);
	void estimatePMD(TParameters & params);
	float calculatePMDS(int readGroup, char & ref, char & read, double & pmdCT, double & pmdGA, double & errorRate, double & pi, float & probPMD, float & probNoPMD);
	void runPMDS(TParameters & params);
	void mergePairedEndReads(TParameters & params);
	void generatePSMCInput(TParameters & params);
	void downSampleBamFile(TParameters & params);
	void estimateApproximateCoverage(TParameters & params);
	void estimateApproximateCoveragePerWindow(TParameters & params);
	void estimateCoveragePerSite(TParameters & params);
	void simulateGWASData(TParameters & params);
	void calculatePoolFreqLikelihoods(TParameters & params);
};




#endif /* LOCI_H_ */
