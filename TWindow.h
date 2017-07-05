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
#include "TReadGroups.h"
#include "bamtools/utils/bamtools_fasta.h"
#include "TRecalibration.h"
#include "TSite.h"
#include "TBedReader.h"
#include "TSiteSubset.h"
#include "TPostMortemDamage.h"

//---------------------------------------------------------------
//EMParameters
//---------------------------------------------------------------
struct EMParameters{
	int numIterations;
	int numThetaOnlyUpdates;
	double maxEpsilon;
	int NewtonRaphsonNumIterations;
	double NewtonRaphsonMaxF;
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
//TWindow
//---------------------------------------------------------------
class TWindow{
public:
	long start;
	long end; //end NOT included in window!
	long length;
	TSite* sites;
	bool sitesInitialized;
	int numReadsInWindow;
	double coverage, fractionSitesNoData, fractionsitesCoverageAtLeastTwo;
	double fractionRefIsN;
	TBaseFrequencies baseFreq;
	TGenotypeMap genoMap;
	bool referenceBaseAdded;

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
	bool addFromRead(BamTools::BamAlignment & bamAlignement, TPMD* pmdObjects, TReadGroups* readGroupsm, int & minQuality, int & maxQuality);
	void addReferenceBaseToSites(BamTools::Fasta & reference, int & refId);
	void addReferenceBaseToSites(TSiteSubset* subset);
	void applyMask(TBedReader* mask, bool inverseMasking);
	void maskCpG(BamTools::Fasta & reference, int & refId);
	void estimateBaseFrequencies();
	void calculateEmissionProbabilities(TRecalibration* recalObject);
	void callMLEGenotype(TRecalibration* recalObject, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF, bool gVCF, bool noAltIfHomoRef);
	void printPileup(TRecalibration* recalObject, std::ofstream & out, std::string & chr);
	void calcCoverage();
	void calcFracN();
	void calcCoveragePerSite(long * siteCoverage, unsigned int maxCov);
	void applyCoverageFilter(int minCoverage, int maxCoverage);
	void compileAlleleicCountsTable(long**** & siteImbalance,long** & siteQuality,  const unsigned int & maxCov);
	double calcLogLikelihood(double* pGenotype);
	void addSitesToBQSR(TRecalibrationBQSR & bqsr, TLog* logfile);
	void addSitesToBQSR(TRecalibrationBQSR & bqsr, TSiteSubset* subset, TLog* logfile);
	void addSitesToQualityTransformTable(TRecalibration* recalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile);
	void addSitesToQualityTransformTable(TRecalibration* recalObject, TRecalibration* otherRecalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile);
	void addSitesToPMDTable(TPMDTables & pmdTables, TLog* logfile);
};

class TWindowDiploid:public TWindow{
protected:
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
	void estimateTheta(EMParameters & constants, TRecalibration* recalObject, std::ofstream & out, TLog* logfile, bool & considerRegions);
	void setTheta(double theta){thetaContainer.setTheta(theta);};
	void calcLikelihoodSurface(TRecalibration* recalObject, std::ofstream & out, int & steps);
	void callMLEGenotypeKnownAlleles(TRecalibration* recalObject, TSiteSubset* subset, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool & isVCF, bool & noAltIfHomoRef, bool & beagle, bool & printOnlyGL);
	void callBayesianGenotype(TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF);
	void callBayesianGenotypeKnownAlleles(TSiteSubset* subset, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr ,bool isVCF);
	void callAllelePresence(TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF, bool noAltIfHomoRef);
	void callAllelePresenceKnwonAlleles(TSiteSubset* subset, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool isVCF, bool noAltIfHomoRef);
	void callRandomBase(TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll);
	void generatePSMCInput(int & blockSize, double & confidence, std::ofstream & out, int & nCharOnLine);
	void addSitesWithDepthTwoOrMoreToVector(std::vector<TSiteDiploid*> & siteVec);
};

//decide which positions to add on the fly
class TWindowDiploidSpecificSites:public TWindowDiploid{
protected:
	long nextId;

public:
	TWindowDiploidSpecificSites(std::vector<TSiteDiploid*> & siteVec);
};

//provide regions file with positions
class TWindowDiploidSiteSubset:public TWindowDiploid{
protected:
	TBedReader* subset;
	long nextId;

public:
	TWindowDiploidSiteSubset(TBedReader* Subset);
	void copySites(TWindowDiploid* other);
};

class TWindowHaploid:public TWindow{
private:
	void fillPGenotype(double* pGenotype);

public:
	TWindowHaploid():TWindow(){};
	TWindowHaploid(long Start, long End):TWindow(Start, End){};
	void initSites(long newLength);
	double calcLogLikelihood();
	void addToRecalibrationEM(TRecalibrationEM & recalObject);
	void addToRecalibrationEM(TRecalibrationEM & recalObject, TSiteSubset* subset);
	void addToExpectedBaseCounts(TRecalibration* recalObject, double** expectedCounts);
	void calculatePoolFreqLikelihoods(int & numChromosomes, Base** majorMinor, gz::ogzstream & out, std::string & chr, bool printAll);
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
	void addToCur(BamTools::BamAlignment & bamAlignement, TPMD* pmdObjects, TReadGroups* readGroups, int & minQuality, int & maxQuality){
		curPointer->addFromRead(bamAlignement, pmdObjects, readGroups, minQuality, maxQuality);
	};
	void addToNext(BamTools::BamAlignment & bamAlignement, TPMD* pmdObjects, TReadGroups* readGroups, int & minQuality, int & maxQuality){
		nextPointer->addFromRead(bamAlignement, pmdObjects, readGroups, minQuality, maxQuality);
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
