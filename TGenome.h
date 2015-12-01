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
	TRecalibration* recalObject;
	bool doRecalibration;
	BamTools::BamReader bamReader;
	BamTools::BamRegion bamRegion;
 	BamTools::SamHeader bamHeader;
 	BamTools::BamAlignment bamAlignement;
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
 	long limitWindows;
 	int limitChr; //for debugging

public:
	TGenome(TLog* Logfile, TParameters & params);
	~TGenome(){
		if(doMasking) delete mask;
		if(fastaReference) reference.Close();
		if(doRecalibration) delete recalObject;
		delete[] pmdObjects;
		if(randomGeneratorInitialized) delete randomGenerator;
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
	void fillSequence(std::vector<double> & vec, std::string & str);
	void calculateLikelihoodSurfaceErrorCalibrationEM(TParameters & params);
	void BQSR(TParameters & params);
	void printQualityTransformation(TParameters & params);
	void recalibrateBamFile(TParameters & params);
	void splitSingleEndReadGroups(TParameters & params);
	void estimatePMD(TParameters & params);
};




#endif /* LOCI_H_ */
