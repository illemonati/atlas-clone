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
//TRecalSearch
//---------------------------------------------------------------
class TRecalSearch{
public:
	double best;
	double min, max;
	double initialMin, initialMax;
	double range;
	double reductionFactor;
	double* rangeSteps;
	int steps, numRangeSteps;
	double* search;
	double* LL;
	bool active;
	bool changed;

	TRecalSearch(double Min, double Max, int Steps);
	TRecalSearch(double Min, double Max, int Steps, double init);
	~TRecalSearch(){
		delete[] search;
		delete[] LL;
		delete[] rangeSteps;
	};
	void initialize(double & Min, double & Max, int & Steps, double Init);

	void fillSearch();
	bool optimizeNextSearch();
	double& at(int index){
		if(active) return search[index];
		else return best;
	};
	void addLL(double value, int index){
		LL[index] += value;
	};
	double& atLL(int index){
		return LL[index];
	};
};

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
	void estimateErrorCalibrationEM(TParameters & params);
	void fillSequence(std::vector<double> & vec, std::string & str);
	void calculateLikelihoodSurfaceErrorCalibrationEM(TParameters & params);
};




#endif /* LOCI_H_ */
