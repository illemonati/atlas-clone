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
#include "TAlignmentParser.h"
#include "TRecalibration.h"
#include "TThetaEstimator.h"
#include "TBedReader.h"
#include "TSiteSubset.h"
#include "TPostMortemDamage.h"
#include "TGLF.h"


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
	long numSitesWithData;
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
	bool addFromRead(TAlignmentParser & alignemntParser, TPMD* pmdObjects, TReadGroups* readGroupsm, int & minQuality, int & maxQuality);
	void addReferenceBaseToSites(BamTools::Fasta & reference, int & refId);
	void addReferenceBaseToSites(TSiteSubset* subset);
	void applyMask(TBedReader* mask, bool inverseMasking);
	void maskCpG(BamTools::Fasta & reference, int & refId);
	void estimateBaseFrequencies();
	void calculateEmissionProbabilities(TRecalibration* recalObject);
	void callMLEGenotype(TRecalibration* recalObject, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF, bool gVCF, bool noAltIfHomoRef);
	void printPileup(TRecalibration* recalObject, std::ofstream & out, std::string & chr);
	virtual void calcCoverage();
	void calcFracN();
	void calcCoveragePerSite(long * siteCoverage, size_t maxCov);
	void applyCoverageFilter(int minCoverage, size_t maxCoverage);
	void createDepthMask(size_t minCoverage, size_t maxCoverage, std::ofstream & outputMaskFile, std::string & chr);
	void addSitesToBQSR(TRecalibrationBQSR & bqsr, TLog* logfile);
	void addSitesToBQSR(TRecalibrationBQSR & bqsr, TSiteSubset* subset, TLog* logfile);
	void addSitesToQualityTransformTable(TRecalibration* recalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile);
	void addSitesToQualityTransformTable(TRecalibration* recalObject, TRecalibration* otherRecalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile);
	void addSitesToPMDTable(TPMDTables & pmdTables, TLog* logfile);
};

class TWindowDiploid:public TWindow{
protected:
	/*
	Theta thetaContainer;

	void fillPGenotype(double* pGenotype, double & expTheta);
	virtual void fillP_G(double* P_g, double* pGenotype);
	virtual double calcLogLikelihood(double* pGenotype);
	void findGoodStartingTheta(Theta & thetaContainer, EMParameters & EMParams);
	void runEMForTheta(Theta & thetaContainer, EMParameters & constants, long & lengthWithData);
	void estimateConfidenceInterval(Theta & thetaContainer);
	*/

public:
	TWindowDiploid():TWindow(){};
	TWindowDiploid(long Start, long End):TWindow(Start, End){};
	void initSites(long newLength);
	void addSitesToThetaEstimator(TRecalibration* recalObject, TThetaEstimator & estimator);
	void addSitesToThetaEstimator(TThetaEstimator & estimator);

	/*
	void estimateTheta(EMParameters & EMParams, TRecalibration* recalObject, std::ofstream & out, TLog* logfile, bool & considerRegions);
	void setTheta(double theta){thetaContainer.setTheta(theta);};
	void calcLikelihoodSurface(TRecalibration* recalObject, std::ofstream & out, int & steps);
	*/

	void callMLEGenotypeKnownAlleles(TRecalibration* recalObject, TSiteSubset* subset, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool & isVCF, bool & noAltIfHomoRef, bool & beagle, bool & printOnlyGL);
	void callBayesianGenotype(TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF);
	void callBayesianGenotypeKnownAlleles(TSiteSubset* subset, TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr ,bool isVCF);
	void callAllelePresence(TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF, bool noAltIfHomoRef);
	void callAllelePresenceKnwonAlleles(TSiteSubset* subset, TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool isVCF, bool noAltIfHomoRef);
	void callRandomBase(TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll);
	void majorityCall(TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll);
	void addToGLF(TGlfWriter & writer, bool printAll);
	void generatePSMCInput(TThetaEstimator & estimator, int & blockSize, double & confidence, std::ofstream & out, int & nCharOnLine);
	void addSitesWithDepthTwoOrMoreToVector(std::vector<TSiteDiploid*> & siteVec);
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
