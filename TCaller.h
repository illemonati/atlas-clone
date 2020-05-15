/*
 * TCaller.h
 *
 *  Created on: Nov 17, 2018
 *      Author: phaentu
 */

#ifndef TCALLER_H_
#define TCALLER_H_

#include "TQualityMap.h"
#include "TRandomGenerator.h"
#include "gzstream.h"
#include "TSite.h"
#include "Vcf/TVCFFields.h"
#include "TLog.h"
#include "TGenotypeData.h"

using namespace GenotypeLikelihoods;

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
	bool _printAltIfHomoRef;
	bool _allowTriallelicSites;
	std::string missingGenotype;

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
	TGenotypeData* genotypePrior; //for callers using a prior. Note: all callers accept priors, but may not use them.
	bool priorSet;

	//functions regarding VCF file
	void setAcceptableFields(TVCFFieldVector* fields, std::string tags);
	void printField(TVCFFieldVector* fields, std::string tag);
	void writeVCFHeader(const std::string & sampleName);

	//function to write info fields
	std::vector< std::string (TCaller::*)(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) > VCFInfoFunctionsVec;
	void fillInfoFieldFunctionPointers();
	virtual std::string getVCFInfoString_DP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

	//functions to write genotype fields
	std::vector< std::string (TCaller::*)(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) > VCFGenotypeFunctionsVec;
	void fillGenotypeFieldFunctionPointers();
	virtual std::string getVCFGenotypeString_GT(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	virtual std::string getVCFGenotypeString_DP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	virtual std::string getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw "Function std::string getVCFGenotypeString_GQ(const const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"; };
	virtual std::string getVCFGenotypeString_AD(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	virtual std::string getVCFGenotypeString_AP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw "Function std::string getVCFGenotypeString_AP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"; };
	virtual std::string getVCFGenotypeString_GL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw "Function std::string getVCFGenotypeString_GL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"; };
	virtual std::string getVCFGenotypeString_PL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw "Function std::string getVCFGenotypeString_PL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"; };
	virtual std::string getVCFGenotypeString_GP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw "Function std::string getVCFGenotypeString_GP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"; };
	virtual std::string getVCFGenotypeString_AB(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw "Function std::string getVCFGenotypeString_AB(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"; };
	virtual std::string getVCFGenotypeString_AI(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw "Function std::string getVCFGenotypeString_AI(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"; };

	//write VCF
	std::string composeVCFString(std::vector<std::string (TCaller::*)(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods)> & vec, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	virtual void writeAlternativeAllelesToVCF();
	void writeCallToVCF(const std::string & chr, const long pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void writeMissingDataToVCF(const TSite & site);
	virtual void clearAfterCall();

	//call
	void countAlleles(const TSite & site);
	virtual void callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	virtual void callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
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
	void setPrintAltIfHomoRef(bool print){ _printAltIfHomoRef = print; };
	void setAllowTriallelic(bool allowTriallelic){ _allowTriallelicSites = allowTriallelic; };
	bool printSitesWithNoData(){ return _printSitesWithNoData; };
	void reportSettings(TLog* logfile);

	//prior
	bool usesPrior(){ return _usesPrior; };
	void setPrior(TGenotypeData* prior){ genotypePrior = prior; priorSet = true; };
	void call(const std::string & chr, const long pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void call(const std::string & chr, const long pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods, const char & firstAllele, const char & secondAllele);
};



//------------------------------------------------------
// TCallerRandomBase
//------------------------------------------------------
class TCallerRandomBase:public TCaller{
private:
	void callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCallerRandomBase(TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerMajorityCall
//------------------------------------------------------
class TCallerMajorityBase:public TCaller{
private:
	void callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCallerMajorityBase(TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerAllelePresence
//------------------------------------------------------
class TCallerAllelePresence:public TCaller{
private:
	TGenotypePosteriorProbabilities posterior;
	double allelePostProb[4];
	int MAP;

	void fillPosteriors(TGenotypeLikelihoods & genotypeLikelihoods);
	void callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string getVCFGenotypeString_AP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

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
	TGenotypeData tmpGenoData;

	void clearAfterCall();
	void callGenotypeFromMetric(TGenotypeData & metric);
	void callGenotypeFromMetricKnownAlleles(const TGenotypeData & metric, std::vector<int> & indeces);
	void callGenotypeFromMetricKnownAllelesUpdateIndex(const TGenotypeData & metric);
	std::string getPerGenotypeMetricString(TGenotypeData & metric);
	void calculateImbalance(const TSite & site);
	std::string getVCFGenotypeString_AB(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string getVCFGenotypeString_AI(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCallerDiploid(TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerMLE
//------------------------------------------------------
class TCallerMLE:public TCallerDiploid{
private:
	void callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string getVCFGenotypeString_GL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string getVCFGenotypeString_PL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCallerMLE(TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerBayes
//------------------------------------------------------
class TCallerBayes:public TCallerDiploid{
private:
	TGenotypePosteriorProbabilities posterior;

	void callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string getVCFGenotypeString_GP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string getVCFGenotypeString_PP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCallerBayes(TRandomGenerator* RandomGenerator);
};

#endif /* TCALLER_H_ */
