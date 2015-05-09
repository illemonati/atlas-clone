/*
 * loci.h
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#ifndef LOCI_H_
#define LOCI_H_

#include "TLog.h"
#include "TParameters.h"
#include "bamtools/api/BamReader.h"
#include "bamtools/api/SamSequenceDictionary.h"

#include "TSite.h"
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>


//---------------------------------------------------------------
//Constants
//---------------------------------------------------------------
struct Constants{
	//post mortem damage
	TPMD pmd_C_T;
	TPMD pmd_G_A;

	//EM parameters
	int numIterations;
	int numThetaOnlyUpdates;
	double maxEpsilon;
	int NewtonRalphsonNumIterations;
	double NewtonRalphsonMaxF;
	double initalTheta;
	double initThetaSearchFactor;
	int initThetaNumSearchIterations;

	//genotype map for enum type
	Genotype** genotypeMap; //mapping base numbering to genotype numbering

	Constants(){
		numIterations = -1;
		numThetaOnlyUpdates = -1;
		maxEpsilon = 0.0;
		NewtonRalphsonNumIterations = -1;
		NewtonRalphsonMaxF = 0.0;
		initalTheta = 0.0;
		initThetaSearchFactor = -1;
		initThetaNumSearchIterations = -1;

		//set up genotype map
		genotypeMap = new Genotype*[4];
		for(int i=0; i<4; ++i)
			genotypeMap[i] = new Genotype[4];

		int geno = 0;
		for(int i=0; i<4; ++i){
			for(int j=i; j<4; ++j){
				genotypeMap[i][j] = static_cast<Genotype>(geno);
				genotypeMap[j][i] = genotypeMap[i][j];
				++geno;
			}
		}
	};

	~Constants(){
		for(int i=0; i<4; ++i) delete[] genotypeMap[i];
		delete[] genotypeMap;
	};

	void initPmd(TParameters & params, TLog* logfile);
	Genotype getGenotype(Base first, Base second){
		return genotypeMap[first][second];
	};
	Genotype getGenotype(int first, int second){
		return genotypeMap[first][second];
	};
	std::string getGenotypeString(int num){
		if(num==0) return "AA";
		if(num==1) return "AC";
		if(num==2) return "AG";
		if(num==3) return "AT";
		if(num==4) return "CC";
		if(num==5) return "CG";
		if(num==6) return "CT";
		if(num==7) return "GG";
		if(num==8) return "GT";
		if(num==9) return "TT";
		return "ERROR!";
	};
};

//---------------------------------------------------------------
//TWindow
//---------------------------------------------------------------
class TWindow{
public:
	long start;
	long end; //end NOT included in window!
	long length;
	TSite* sites;
	bool sitesInitialized;
	double coverage, fractionSitesNoData, fractionsitesCoverageAtLeastTwo;;
	TBaseFrequencies baseFreq;
	double theta;
	double LL;
	Constants* sharedConstants;
	TWindow(Constants* SharedConstants);
	TWindow(long Start, long End, Constants* SharedConstants);
	~TWindow(){
		if(sitesInitialized) delete[] sites;
	};
	void initSites(long newLength);
	void clear();
	void move(long Start, long End);
	bool addFromRead(BamTools::BamAlignment & bamAlignement);
	void fillPGenotype(double* pGenotype, double & expTheta);
	void fillP_G(double* P_g, double* pGenotype);
	void calcLogLikelihood(double* pGenotype);
	void findGoodStartingTheta();
	void runEM(TLog* logfile);
	void writeEMResults(std::ofstream & out);
	void printPileup(std::ofstream & out, std::string & chr);
	void calcCoverage();
	void calcLikelihoodSurface(std::ofstream & out, std::string & chr);
};

//---------------------------------------------------------------
//TGenome
//---------------------------------------------------------------
class TGenome{
public:
	//std::vector<TChromosome> chromosomes;
	std::string filename;
	std::string outname;
	TLog* logfile;
	BamTools::BamReader bamReader;
 	BamTools::BamRegion bamRegion;
 	BamTools::SamHeader bamHeader;
 	BamTools::BamAlignment bamAlignement;
 	BamTools::SamSequenceIterator chrIterator;
 	int chrNumber;
 	long chrLength;
 	long curStart, curEnd;
 	long windowSize;
 	double maxMissing;
 	Constants sharedConstants;
 	long oldPos;
 	TWindow* curWindow;
 	TWindow* nextWindow;
 	std::string outputName;

	TGenome(TLog* Logfile, TParameters & params);
	~TGenome(){
		delete curWindow;
		delete nextWindow;
	};
	bool iterateChromosome();
	bool iterateWindow();
	bool readData();
	void estimateTheta();
	void calcLikelihoodSurfaces();
	void printPileup();
};




#endif /* LOCI_H_ */
