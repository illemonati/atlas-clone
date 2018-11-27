/*
 * TCaller.h
 *
 *  Created on: Nov 17, 2018
 *      Author: phaentu
 */

#ifndef TCALLER_H_
#define TCALLER_H_

#include "TGenotypeMap.h"
#include "TRandomGenerator.h"
#include "gzstream.h"
#include "TSite.h"
#include "TVCFFields.h"
#include "TLog.h"

//------------------------------------------------------
// TCaller
// Note: this is base class, not meant to be used but to derive from
//------------------------------------------------------
class TCaller{
protected:
	//caller specific defaults
	std::string callerName;
	std::string filenameExtention;

	//lookup stuff
	TGenotypeMap genoMap;
	TQualityMap qualMap;
	TVCFInfoFields VCFInfoFields;
	TVCFGenotypeFields VCFGenotypeFields;
	TRandomGenerator* randomGenerator;

	//output choices
	bool _printSitesWithNoData;
	bool _printNoAltIfHomoRef;

	//output file
	std::string filename;
	gz::ogzstream vcf;
	bool vcfOpen;
	std::string genotypeFormatString;

	//temp variables for calling
	std::string calledGenotype;
	std::vector<int> genotypesWithHighestMetric;
	int referenceBase;
	std::vector<int> altAlleles; //order of Base enums: A, C, G, T, N
	int alleleCounts[4];
	bool allelesCounted;

	//genotype prior
	bool _usesPrior;
	double* genotypePrior; //for callers using a prior. Note: all callers accept priors, but may not use them.
	bool priorSet;

	//functions regarding VCF file
	void setAcceptableFields(TVCFFieldVector* fields, std::string tags);
	void printField(TVCFFieldVector* fields, std::string tag);
	void writeVCFHeader(const std::string & sampleName);

	//function to write info fields
	std::vector<std::string (TCaller::*)(TSite & site)> VCFInfoFunctionsVec;
	void fillInfoFieldFunctionPointers();
	virtual std::string getVCFInfoString_DP(TSite & site);

	//functions to write genotype fields
	std::vector<std::string (TCaller::*)(TSite & site)> VCFGenotypeFunctionsVec;
	void fillGenotypeFieldFunctionPointers();
	virtual std::string getVCFGenotypeString_GT(TSite & site);
	virtual std::string getVCFGenotypeString_DP(TSite & site);
	virtual std::string getVCFGenotypeString_GQ(TSite & site){ throw "Function std::string getVCFGenotypeString_GQ(const TSite & site) not defined for base class TCaller!"; };
	virtual std::string getVCFGenotypeString_AD(TSite & site);
	virtual std::string getVCFGenotypeString_AP(TSite & site){ throw "Function std::string getVCFGenotypeString_AP(const TSite & site) not defined for base class TCaller!"; };
	virtual std::string getVCFGenotypeString_GL(TSite & site){ throw "Function std::string getVCFGenotypeString_GL(const TSite & site) not defined for base class TCaller!"; };
	virtual std::string getVCFGenotypeString_PL(TSite & site){ throw "Function std::string getVCFGenotypeString_PL(const TSite & site) not defined for base class TCaller!"; };
	virtual std::string getVCFGenotypeString_GP(TSite & site){ throw "Function std::string getVCFGenotypeString_GP(const TSite & site) not defined for base class TCaller!"; };
	virtual std::string getVCFGenotypeString_AB(TSite & site){ throw "Function std::string getVCFGenotypeString_AB(const TSite & site) not defined for base class TCaller!"; };
	virtual std::string getVCFGenotypeString_AI(TSite & site){ throw "Function std::string getVCFGenotypeString_AI(const TSite & site) not defined for base class TCaller!"; };

	//write VCF
	std::string composeVCFString(std::vector<std::string (TCaller::*)(TSite & site)> & vec, TSite & site);
	virtual void writeAlternativeAllelesToVCF();
	void writeCallToVCF(const std::string & chr, const long pos, TSite & site);
	virtual void clearAfterCall();

	//call
	void countAlleles(TSite & site);
	virtual void callGenotype(TSite & site);
	virtual void callGenotypeKnownAlleles(TSite & site);
	template <typename T> int pickIndexWithHighestMetric(T* metric, const int size);
	template <typename T> int pickIndexWithSecondHighestMetric(T* metric, const int size, const int excludeIndex);

public:
	TCaller(TRandomGenerator* RandomGenerator);
	virtual ~TCaller();

	//set which fields to print
	void printInfoFields(std::vector<std::string> & tags);
	void printInfoFields(std::string tags);
	void printGenotypeFields(std::vector<std::string> & tags);
	void printGenotypeFields(std::string tags);

	//open / close vcf file
	void openVCF(const std::string Filename, const std::string sampleName);
	void closeVCF();

	//other output options
	void setPrintSitesWithNoData(bool print){ _printSitesWithNoData = print; };
	void setNoAltIfHomoRef(bool print){ _printNoAltIfHomoRef = print; };
	bool printSitesWithNoData(){ return _printSitesWithNoData; };
	void reportSettings(TLog* logfile);

	//prior
	bool usesPrior(){ return _usesPrior; };
	void setPrior(double* GenoPrior){ genotypePrior = GenoPrior; priorSet = true; };
	void call(const std::string & chr, const long pos, TSite & site);
	void call(const std::string & chr, const long pos, TSite & site, char & first, char & second);
};



//------------------------------------------------------
// TCallerRandomBase
//------------------------------------------------------
class TCallerRandomBase:public TCaller{
private:
	void callGenotype(TSite & site);
	void callGenotypeKnownAlleles(TSite & site);

public:
	TCallerRandomBase(TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerMajorityCall
//------------------------------------------------------
class TCallerMajorityBase:public TCaller{
private:
	void callGenotype(TSite & site);
	void callGenotypeKnownAlleles(TSite & site);

public:
	TCallerMajorityBase(TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerAllelePresence
//------------------------------------------------------
class TCallerAllelePresence:public TCaller{
private:
	double posteriorProb[10];
	double allelePostProb[4];
	int MAP;

	void callGenotype(TSite & site);
	std::string getVCFGenotypeString_GQ(TSite & site);
	std::string getVCFGenotypeString_AP(TSite & site);

public:
	TCallerAllelePresence(TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerDiploid
//------------------------------------------------------
class TCallerDiploid:public TCaller{
protected:
	int indexOfMax, indexOfSecond;
	std::string AB, AI;
	bool imbalanceCalculated;

	void clearAfterCall();
	void callGenotypeFromMetric(double* metric);
	void callGenotypeFromMetricKnownAlleles(double* metric);
	std::string getPerGenotypeMetricString(double* metric);
	void calculateImbalance(TSite & site);
	std::string getVCFGenotypeString_AB(TSite & site);
	std::string getVCFGenotypeString_AI(TSite & site);

public:
	TCallerDiploid(TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerMLE
//------------------------------------------------------
class TCallerMLE:public TCallerDiploid{
private:
	void callGenotype(TSite & site);
	void callGenotypeKnownAlleles(TSite & site);
	std::string getVCFGenotypeString_GQ(TSite & site);
	std::string getVCFGenotypeString_GL(TSite & site);
	std::string getVCFGenotypeString_PL(TSite & site);

public:
	TCallerMLE(TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerBayes
//------------------------------------------------------
class TCallerBayes:public TCallerDiploid{
private:
	double posteriorProb[10];

	void callGenotype(TSite & site);
	void callGenotypeKnownAlleles(TSite & site);
	std::string getVCFGenotypeString_GQ(TSite & site);
	std::string getVCFGenotypeString_GP(TSite & site);
	std::string getVCFGenotypeString_PP(TSite & site);

public:
	TCallerBayes(TRandomGenerator* RandomGenerator);
};

#endif /* TCALLER_H_ */
