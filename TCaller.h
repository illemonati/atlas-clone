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

	//std::vector<std::string> acceptableVCFInfoFields;
	//std::vector<VCFGenotypeFieldTag> acceptableVCFGenotypeFields;

	//lookup stuff
	TGenotypeMap genoMap;
	TQualityMap qualMap;
	TVCFInfoFields VCFInfoFields;
	TVCFGenotypeFields VCFGenotypeFields;
	TRandomGenerator* randomGenerator;

	//output choices
	bool _printSitesWithNoData;

	//output file
	std::string filename;
	gz::ogzstream vcf;
	bool vcfOpen;
	std::string genotypeFormatString;

	//temp variables for calling
	int numGenotypes;
	double perGenotypeMetric[10]; //could be likelihood or posterior prob per genotype, depending on metric
	std::string calledGenotype;
	std::vector<int> genotypesWithHighestMetric;
	std::vector<char> altAlleles;

	//functions regarding VCF file
	void setAcceptableFields(TVCFFieldVector* fields, std::string tags);
	void printField(TVCFFieldVector* fields, std::string tag);
	void writeVCFHeader(const std::string & sampleName);

	//calling
	virtual void fillCallingMetric(TSite & site){ throw "fillCallingMetric() not defined for base class TCaller!"; };
	virtual void callGenotype(TSite & site);


	//function to write info fields
	std::vector<std::string (TCaller::*)(const TSite & site)> VCFInfoFunctionsVec;
	void fillInfoFieldFunctionPointers();
	virtual std::string getVCFInfoString_DP(const TSite & site);

	//functions to write genotype fields
	std::vector<std::string (TCaller::*)(const TSite & site)> VCFGenotypeFunctionsVec;
	void fillGenotypeFieldFunctionPointers();
	virtual std::string getVCFGenotypeString_GT(const TSite & site);
	virtual std::string getVCFGenotypeString_DP(const TSite & site);
	virtual std::string getVCFGenotypeString_GQ(const TSite & site){ throw "Function std::string getVCFGenotypeString_GQ(const TSite & site) not defined for base class TCaller!"; };

	//write VCF
	std::string composeVCFString(std::vector<std::string (TCaller::*)(const TSite & site)> & vec, const TSite & site);
	void writeCallToVCF(const std::string & chr, const long pos, TSite & site);


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
// TCallerDiploid
//------------------------------------------------------
class TCallerDiploid:public TCaller{
protected:
	//output choices
	bool peformImbalanceTest;

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
class TCallerBayes:public TCaller{
private:
	double genotypePriorProbabilities[10];

	void fillCallingMetric(TSite & site);

public:
	TCallerBayes(TRandomGenerator* RandomGenerator);

	void setPrior(double* GenotypePriorProbabilities);

};

#endif /* TCALLER_H_ */
