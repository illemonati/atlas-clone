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

//---------------------------------------------------------------
//TGenome
//---------------------------------------------------------------
class TGenome{
private:
	void restartChromosome(TWindowPair & windowPair);
	bool iterateChromosome(TWindowPair & windowPair);
	void moveChromosome(TWindowPair & windowPair);
	bool iterateWindow(TWindowPair & windowPair);
	bool readData(TWindowPair & windowPair);
	void initializePostMortemDamage(TParameters & params);
	void initializeRecalibration(TParameters & params);
	void openThetaOutputFile(std::ofstream & out);

	TPMD pmdObject;
	TRecalibration* recalObject;
	BamTools::BamReader bamReader;
	BamTools::BamRegion bamRegion;
 	BamTools::SamHeader bamHeader;
 	BamTools::BamAlignment bamAlignement;
 	BamTools::SamSequenceIterator chrIterator;
 	TReadGroups readGroups;
 	int chrNumber;
 	long chrLength;
 	long curStart;
 	long curEnd;
	std::string filename;
	TLog* logfile;
 	int windowSize;
 	int numWindowsOnChr;
 	double maxMissing;
 	long oldPos;
 	std::string outputName;
 	TBedReader* mask;
 	bool doMasking;

public:
	TGenome(TLog* Logfile, TParameters & params);
	~TGenome(){
		if(doMasking) delete mask;
		delete recalObject;
	};
	void estimateTheta(TParameters & params);
	void calcLikelihoodSurfaces(TParameters & params);
	void callMLEGenotypes(TParameters & params);
	void callAllelePresence(TParameters & params);
	void printPileup();
	void estimateErrorCalibration(TParameters & params);
	void estimateErrorCalibrationEM(TParameters & params);
	void fillSequence(std::vector<double> & vec, std::string & str);
	void calculateLikelihoodSurfaceErrorCalibrationEM(TParameters & params);
	void BQSR(TParameters & params);
};




#endif /* LOCI_H_ */
