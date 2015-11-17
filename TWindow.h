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
#include "bamtools/utils/bamtools_fasta.h"
#include "TRecalibration.h"
#include "TSite.h"
#include "TBedReader.h"
#include "TSiteSubset.h"

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
	EMParameters(TParameters & params, TLog* logfile);
	~EMParameters(){};

	void report(TLog* logfile);
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
//TReadGroups
//---------------------------------------------------------------
struct readGroup{
public:
	std::string name;
	int id;
	BamTools::SamReadGroup* object;
};

class TReadGroups{
public:
	readGroup* groups;
	int numGroups;
	bool initialized;

	TReadGroups(){
		initialized = false;
		numGroups = 0;
		groups = NULL;
	};

	~TReadGroups(){
		if(initialized) delete[] groups;
	};

	void fill(BamTools::SamHeader & bamHeader){
		//empty if filled before
		if(initialized) delete[] groups;
		//create and fill array
		numGroups = bamHeader.ReadGroups.Size();
		groups = new readGroup[numGroups];
		int i = 0;
		for(BamTools::SamReadGroupIterator it = bamHeader.ReadGroups.Begin(); it != bamHeader.ReadGroups.End(); ++it, ++i){
			groups[i].id = i;
			groups[i].name = it->ID;
			groups[i].object= &(*it);
		}
		initialized = true;
	};

	int find(std::string & name){
		for(int i=0; i<numGroups; ++i){
			if(groups[i].name == name) return i;
		}
		throw "Read Group '" + name + "' was not present in header of bam file!";
	};

	bool readGroupExists(std::string & name){
		for(int i=0; i<numGroups; ++i){
			if(groups[i].name == name) return true;
		}
		return false;
	};

	std::string getName(int num){
		if(num < 0 || num >= numGroups) throw "No read group with number " + toString(num) + "!";
		return groups[num].name;
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
	double coverage, fractionSitesNoData, fractionsitesCoverageAtLeastTwo;
	TBaseFrequencies baseFreq;
	TGenotypeMap genoMap;

	TWindow();
	TWindow(long Start, long End);
	virtual ~TWindow(){
		if(sitesInitialized) delete[] sites;
	};
	virtual void initSites(long newLength){
		throw "Function 'initSites' not implemented for base class TWindow!";
	};
	void clear();
	void move(long Start, long End);
	bool addFromRead(BamTools::BamAlignment & bamAlignement, TPMD* pmdObjects, TReadGroups* readGroups);
	void addReferenceBaseToSites(BamTools::Fasta & reference, int & refId);
	void addReferenceBaseToSites(TSiteSubset* subset);
	void applyMask(TBedReader* mask);
	void estimateBaseFrequencies();
	void calculateEmissionProbabilities(TRecalibration* recalObject);
	void callMLEGenotype(TRecalibration* recalObject, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF);
	void callMLEGenotypeKnownAlleles(TRecalibration* recalObject, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr);
	void printPileup(TRecalibration* recalObject, std::ofstream & out, std::string & chr);
	void calcCoverage();
	double calcLogLikelihood(double* pGenotype);
	void addSitesToBQSR(TRecalibrationBQSR & bqsr, TLog* logfile);
	void addSitesToQualityTransformTable(TRecalibration* recalObject, TQualityTransformTable & QT, TLog* logfile);
};

class TWindowDiploid:public TWindow{
private:
	Theta thetaContainer;

	void fillPGenotype(double* pGenotype, double & expTheta);
	void fillP_G(double* P_g, double* pGenotype);
	void findGoodStartingTheta(Theta & thetaContainer, EMParameters & EMParams);
	void runEMForTheta(Theta & thetaContainer, EMParameters & constants, int & lengthWithData);
	void estimateConfidenceInterval(Theta & thetaContainer);

public:
	TWindowDiploid():TWindow(){};
	TWindowDiploid(long Start, long End):TWindow(Start, End){};
	void initSites(long newLength);
	void estimateTheta(EMParameters & constants, TRecalibration* recalObject, std::ofstream & out, TLog* logfile);
	void setTheta(double theta){thetaContainer.setTheta(theta);};
	void calcLikelihoodSurface(TRecalibration* recalObject, std::ofstream & out, int & steps);
	void callAllelePresence(TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF);
	void callBayesianGenotype(TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF);
	void callBayesianGenotypeKnownAlleles(TSiteSubset* subset, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printRef, bool isVCF);
};

class TWindowHaploid:public TWindow{
private:
	void fillPGenotype(double* pGenotype);

public:
	TWindowHaploid():TWindow(){};
	TWindowHaploid(long Start, long End):TWindow(Start, End){};
	void initSites(long newLength);
	double calcLogLikelihood();
	void addToJacobian(TRecalibrationEM* reclObject);
	void addToLikelihoodRecalibration(TRecalibrationEM* reclObject);
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
	void addToCur(BamTools::BamAlignment & bamAlignement, TPMD* pmdObjects, TReadGroups* readGroups){
		curPointer->addFromRead(bamAlignement, pmdObjects, readGroups);
	};
	void addToNext(BamTools::BamAlignment & bamAlignement, TPMD* pmdObjects, TReadGroups* readGroups){
		nextPointer->addFromRead(bamAlignement, pmdObjects, readGroups);
	};
	void clear(){
		curPointer->clear();
		nextPointer->clear();
	}
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
