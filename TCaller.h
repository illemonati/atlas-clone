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
	std::string defaultInfoFields, defaultGenotypeFields;

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
	double genotypePrior[10]; //for callers using a prior. Note: all callers accept priors, but may not use them.
	void setPrior(double* GenotypePrior){ for(int g=0; g<10; ++g) genotypePrior[g] = GenotypePrior[g]; };

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

	//write VCF
	std::string composeVCFString(std::vector<std::string (TCaller::*)(TSite & site)> & vec, TSite & site);
	virtual void writeAlternativeAllelesToVCF();
	void writeCallToVCF(const std::string & chr, const long pos, TSite & site);

	//call
	void countAlleles(TSite & site);
	virtual void callGenotype(TSite & site);
	template <typename T> int pickIndexWithHighestMetric(T* metric, const int size, double & maxMetric);
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


	void call(const std::string & chr, const long pos, TSite & site);
};



//------------------------------------------------------
// TCallerRandomBase
//------------------------------------------------------
class TCallerRandomBase:public TCaller{
private:
	virtual void callGenotype(TSite & site);

public:
	TCallerRandomBase(TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerMajorityCall
//------------------------------------------------------
class TCallerMajorityBase:public TCaller{
private:
	double highestCounts;
	virtual void callGenotype(TSite & site);

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
	double highestPostProb;

	void callGenotype(TSite & site);
	std::string getVCFGenotypeString_GQ(TSite & site);

public:
	TCallerAllelePresence(TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerDiploid
//------------------------------------------------------
class TCallerDiploid:public TCaller{
protected:
	//output choices
	bool peformImbalanceTest;
	int numGenotypes;
	double perGenotypeMetric[10]; //either likelihood or posterior probability

public:
	TCallerDiploid(TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerMLE
//------------------------------------------------------
class TCallerMLE:public TCallerDiploid{
private:
	bool _print_GL;
	bool _print_PL;

	void fillCallingMetric(TSite & site);

public:
	TCallerMLE(TRandomGenerator* RandomGenerator);

	void setPrint_GL(bool print){ _print_GL = print; };
	bool print_GL(){ return _print_GL; };

};

//------------------------------------------------------
// TCallerBayes
//------------------------------------------------------
class TCallerBayes:public TCallerDiploid{
private:
	double genotypePriorProbabilities[10];

	void fillCallingMetric(TSite & site);

public:
	TCallerBayes(TRandomGenerator* RandomGenerator);

	void setPrior(double* GenotypePriorProbabilities);

};

#endif /* TCALLER_H_ */
