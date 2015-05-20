/*
 * TWindow.h
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#ifndef TWINDOW_H_
#define TWINDOW_H_

#include "TLog.h"
#include "TParameters.h"
#include "bamtools/api/BamReader.h"
#include "bamtools/api/SamSequenceDictionary.h"

#include "TSite.h"
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

//---------------------------------------------------------------
//EMParameters
//---------------------------------------------------------------
struct EMParameters{
	int numIterations;
	int numThetaOnlyUpdates;
	double maxEpsilon;
	int NewtonRalphsonNumIterations;
	double NewtonRalphsonMaxF;
	double initalTheta;
	double initThetaSearchFactor;
	int initThetaNumSearchIterations;

	EMParameters();
	~EMParameters(){};
};

//---------------------------------------------------------------
//Recalibration
//---------------------------------------------------------------
struct Recalibration{
	bool doRecalibration;
	double a,b;

	Recalibration(){
		doRecalibration = false;
		a = 0.0;
		b = 0.0;
	};
	Recalibration(const double & paramA, const double & paramB){
		doRecalibration = true;
		a = paramA;
		b = paramB;
	};
	Recalibration(std::string recalString);
	~Recalibration(){};
};

//---------------------------------------------------------------
//Theta
//---------------------------------------------------------------
struct Theta{
	double theta, expTheta, thetaConfidence, LL;

	Theta(){
		theta = 0.0;
		thetaConfidence = 0.0;
		expTheta = 0.0;
		LL = -9e100;
	};

	void setTheta(double val){
		theta = val;
		expTheta = exp(-theta);
	}
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
	double coverage, fractionSitesNoData, fractionsitesCoverageAtLeastTwo;
	TBaseFrequencies baseFreq;


	TWindow();
	TWindow(long Start, long End);
	virtual ~TWindow(){
		if(sitesInitialized) delete[] sites;
	};
	virtual void initSites(long newLength){ throw "Function 'initSites' not implemented for base class TWindow!"; };
	void clear();
	void move(long Start, long End);
	bool addFromRead(BamTools::BamAlignment & bamAlignement);
	void estimateBaseFrequencies();
	void calculateEmissionProbabilities(TPMD & pmdObject);
	void printPileup(TPMD & pmd, std::ofstream & out, std::string & chr);
	void calcCoverage();
	double calcLogLikelihood(double* pGenotype);
	void calculateEissionProbabilitiesWithScaledError(TPMD & pmdObject, double & a, double & b);
};

class TWindowDiploid:public TWindow{
private:
	GenotypeMap genoMap;
	void fillPGenotype(double* pGenotype, double & expTheta);
	void fillP_G(double* P_g, double* pGenotype);
	void findGoodStartingTheta(Theta & thetaContainer, EMParameters & EMParams);
	void runEMForTheta(Theta & thetaContainer, EMParameters & constants, int & lengthWithData);
	void estimateConfidenceInterval(Theta & thetaContainer);

public:
	TWindowDiploid():TWindow(){};
	TWindowDiploid(long Start, long End):TWindow(Start, End){};
	void initSites(long newLength);
	void estimateTheta(EMParameters & constants, TPMD & pmd, Recalibration & recal, std::ofstream & out, TLog* logfile);
	void calcLikelihoodSurface(TPMD & pmd, Recalibration & recal, std::ofstream & out, int & steps);
};

class TWindowHaploid:public TWindow{
private:
	void fillPGenotype(double* pGenotype);

public:
	TWindowHaploid():TWindow(){};
	TWindowHaploid(long Start, long End):TWindow(Start, End){};
	void initSites(long newLength);
	double calcLogLikelihood();
};

//---------------------------------------------------------------
//TWindowPair
//---------------------------------------------------------------
class TWindowPair{
public:
	TWindow* curPointer;
	TWindow* nextPointer;

	TWindowPair(){
		curPointer = NULL;
		nextPointer = NULL;
	};
	virtual ~TWindowPair(){};
	virtual void swap(){
		TWindow* tmp = curPointer;
		curPointer = nextPointer;
		nextPointer = tmp;
	};
	void addToCur(BamTools::BamAlignment & bamAlignement){
		curPointer->addFromRead(bamAlignement);
	};
	void addToNext(BamTools::BamAlignment & bamAlignement){
		nextPointer->addFromRead(bamAlignement);
	};
};

class TWindowPairDiploid:public TWindowPair{
public:
	TWindowDiploid* cur;
	TWindowDiploid* next;

	TWindowPairDiploid(){
		cur = new TWindowDiploid();
		curPointer = cur;
		next = new TWindowDiploid();
		nextPointer = next;
	};
	~TWindowPairDiploid(){
		delete cur;
		delete next;
	};

	void swap(){
		TWindowDiploid* tmp = cur;
		cur = next;
		next = tmp;
		TWindowPair::swap();
	};
};

class TWindowPairHaploid:public TWindowPair{
public:
	TWindowHaploid* cur;
	TWindowHaploid* next;

	TWindowPairHaploid(){
		cur = new TWindowHaploid();
		curPointer = cur;
		next = new TWindowHaploid();
		nextPointer = next;
	};
	~TWindowPairHaploid(){
		delete cur;
		delete next;
	};

	void swap(){
		TWindowHaploid* tmp = cur;
		cur = next;
		next = tmp;
		TWindowPair::swap();
	};
};

#endif /* TWINDOW_H_ */
