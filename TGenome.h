/*
 * loci.h
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#ifndef LOCI_H_
#define LOCI_H_

#include "TWindow.h"

//---------------------------------------------------------------
//TGenome
//---------------------------------------------------------------
class TGenome{
private:
	bool iterateChromosome(TWindowPair & windowPair);
	bool iterateWindow(TWindowPair & windowPair);
	bool readData(TWindowPair & windowPair);
	void initializePostMortemDamage(TParameters & params, TLog* logfile);

	TPMD pmdObject;
	BamTools::BamReader bamReader;
	BamTools::BamRegion bamRegion;
 	BamTools::SamHeader bamHeader;
 	BamTools::BamAlignment bamAlignement;
 	BamTools::SamSequenceIterator chrIterator;
 	int chrNumber;
 	long chrLength;
 	long curStart, curEnd;
	std::string filename;
	std::string outname;
	TLog* logfile;
 	long windowSize;
 	double maxMissing;
 	long oldPos;
 	std::string outputName;

public:
	TGenome(TLog* Logfile, TParameters & params);
	~TGenome(){};
	void estimateTheta(TParameters & params);
	void calcLikelihoodSurfaces(TParameters & params);
	void printPileup();
	void estimateErrorCalibration(TParameters & params);
};




#endif /* LOCI_H_ */
