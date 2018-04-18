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
private:
	void fillPGenotype(double* pGenotype);

public:
	long start;
	long end; //end NOT included in window!
	long length;
	TSite* sites;
	bool sitesInitialized;
	int numReadsInWindow;
	double depth, fractionSitesNoData, fractionsitesDepthAtLeastTwo;
	double fractionRefIsN;
	long numSitesWithData;
	TBaseFrequencies baseFreq;
	TGenotypeMap genoMap;
	bool referenceBaseAdded;

	TWindow();
	TWindow(long Start, long End);
	~TWindow(){
		clear();
		if(sitesInitialized){
			delete[] sites;
		}
	};
	void initSites(long newLength);
	void clear();
	void move(long Start, long End);
	bool addFromRead(TAlignmentParser & alignemntParser, TPMD* pmdObjects);
	void addReferenceBaseToSites(BamTools::Fasta & reference, int & refId);
	void addReferenceBaseToSites(TSiteSubset* subset);
	void applyMask(TBedReader* mask, bool inverseMasking);
	void maskCpG(BamTools::Fasta & reference, int & refId);
	void estimateBaseFrequencies();
	void calculateEmissionProbabilities(TRecalibration* recalObject);
	void callMLEGenotype(TRecalibration* recalObject, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF, bool gVCF, bool noAltIfHomoRef);
	void printPileup(TRecalibration* recalObject, gz::ogzstream & out, std::string & chr);
	void calcDepth();
	void calcFracN();
	void calcDepthPerSite(long * siteDepth, size_t maxCov);
	void printDepthPerSite(gz::ogzstream & out, std::string & chr);
	void countAlleles(long**** siteImbalance, const unsigned int & maxDepth);
	void applyDepthFilter(int minDepth, size_t maxDepth);
	void createDepthMask(size_t minDepth, size_t maxDepth, std::ofstream & outputMaskFile, std::string & chr);

	//add sites to data structures
	void addSitesToBQSR(TRecalibrationBQSR & bqsr, TLog* logfile);
	void addSitesToBQSR(TRecalibrationBQSR & bqsr, TSiteSubset* subset, TLog* logfile);
	void addSitesToQualityTransformTable(TRecalibration* recalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile);
	void addSitesToQualityTransformTable(TRecalibration* recalObject, TRecalibration* otherRecalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile);
	void addSitesToPMDTable(TPMDTables & pmdTables, TLog* logfile);
	void addSitesToThetaEstimator(TRecalibration* recalObject, TThetaEstimator & estimator);
	void addSitesToThetaEstimator(TThetaEstimator & estimator);
	void addToGLF(TGlfWriter & writer, bool printAll);
	void addToRecalibrationEM(TRecalibrationEM & recalObject);
	void addToRecalibrationEM(TRecalibrationEM & recalObject, TSiteSubset* subset);


	//callers
	void callMLEGenotypeKnownAlleles(TRecalibration* recalObject, TSiteSubset* subset, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool & isVCF, bool & noAltIfHomoRef, bool & beagle, bool & printOnlyGL);
	void callBayesianGenotype(TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF);
	void callBayesianGenotypeKnownAlleles(TSiteSubset* subset, TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr ,bool isVCF);
	void callAllelePresence(TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF, bool noAltIfHomoRef);
	void callAllelePresenceKnwonAlleles(TSiteSubset* subset, TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool isVCF, bool noAltIfHomoRef);
	void callRandomBase(TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll);
	void majorityCall(TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll);

	//other
	void generatePSMCInput(TThetaEstimator & estimator, int & blockSize, double & confidence, std::ofstream & out, int & nCharOnLine);
	double calcLogLikelihood();
};

//---------------------------------------------------------------
//TWindowPair
//---------------------------------------------------------------


class TWindowPair{
public:
	TWindow* cur;
	TWindow* next;

	TWindowPair(){
		cur = new TWindow();
		next = new TWindow();
	};
	~TWindowPair(){
		delete cur;
		delete next;
	};

	void swap(){
		TWindow* tmp = cur;
		cur = next;
		next = tmp;
	}
	bool addToCur(TAlignmentParser & alignemntParser, TPMD* pmdObjects){
		return cur->addFromRead(alignemntParser, pmdObjects);
	};
	bool addToNext(TAlignmentParser & alignemntParser, TPMD* pmdObjects){
		return next->addFromRead(alignemntParser, pmdObjects);
	};
	void clear(){
		cur->clear();
		next->clear();
	}

};


#endif /* TWINDOW_H_ */
