/*
 * TCaller.h
 *
 *  Created on: Nov 17, 2018
 *      Author: phaentu
 */

#ifndef TCALLER_H_
#define TCALLER_H_

#include "gzstream.h"
#include "TGenotypeData.h"
#include "TGenome.h"
#include "TTask.h"
#include "TVCFFields.h"

namespace GenomeTasks{

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
	VCF::TVCFInfoFields VCFInfoFields;
	VCF::TVCFGenotypeFields VCFGenotypeFields;
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
	std::vector<Base> genotypesWithHighestMetric;
	Base referenceBase;
	std::vector<Base> altAlleles; //order of Base enums: A, C, G, T, N
	TBaseCounts alleleCounts;
	bool allelesCounted;

	//genotype prior
	bool _usesPrior;
	TGenotypeData* genotypePrior; //for callers using a prior. Note: all callers accept priors, but may not use them.
	bool priorSet;

	//functions regarding VCF file
	void setAcceptableFields(VCF::TVCFFieldVector* fields, std::string tags);
	void printField(VCF::TVCFFieldVector* fields, std::string tag);
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

	template <typename T> uint8_t pickIndexWithHighestMetric(T* metric, const uint8_t size){
		//find maximum
		double maxMetric = 0.0;
		for(uint8_t i=0; i<size; ++i){
			if(metric[i] > maxMetric)
				maxMetric = metric[i];
		}

		//get vec of all index at maximum
		std::vector<uint8_t> vec;
		for(uint8_t i=0; i<size; ++i){
			if(metric[i] == maxMetric)
				vec.push_back(i);
		}

		//return random index among those at max
		return vec[randomGenerator->pickOne(vec.size())];
	};

	template <typename T> uint8_t pickIndexWithSecondHighestMetric(T* metric, const uint8_t size, const uint8_t excludeIndex){
		//find maximum
		double max = 0.0;
		for(uint8_t i=0; i<size; ++i){
			if(i != excludeIndex && metric[i] > max)
				max = metric[i];
		}

		//get vec of all index at maximum
		std::vector<uint8_t> vec;
		for(uint8_t i=0; i<size; ++i){
			if(i != excludeIndex && metric[i] == max)
				vec.push_back(i);
		}

		//return random index among those at max
		return vec[randomGenerator->pickOne(vec.size())];
	};

public:
	TCaller(TRandomGenerator* RandomGenerator);
	virtual ~TCaller();

	std::string name() const{ return callerName; };

	//set which fields to print
	void printInfoFields(std::vector<std::string> & tags);
	void printInfoFields(std::string tags);
	void printGenotypeFields(std::vector<std::string> & tags);
	void printGenotypeFields(std::string tags);

	//open / close vcf file
	void openVCF(const std::string Filename, const std::string sampleName, TLog* logfile);
	void closeVCF();

	//other output options
	void initializeOutput(TParameters & Parameters, TLog* Logfile);
	bool printSitesWithNoData(){ return _printSitesWithNoData; };

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
	TGenotypeProbabilities posterior;
	double allelePostProb[4];
	Base MAP;

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
	TBinomPValue _binomP;

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
	TGenotypeProbabilities posterior;

	void callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string getVCFGenotypeString_GP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string getVCFGenotypeString_PP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCallerBayes(TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCall
// the class to perform calls based on windows
//------------------------------------------------------
class TCall:public TGenome_windows{
private:
	TCaller* _caller;
	TGenotypePrior* _prior;

	void _initializeGenotypePrior(TParameters & Parameters);
	void _call();
	void _callKnwonAlleles();
	void _handleWindow();

public:
	TCall(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	~TCall();
	void call();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_call:public TTask{
public:
	TTask_call(){ _explanation = "Calling genotypes"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TCall caller(Parameters, Logfile, _randomGenerator);
		caller.call();
	};
};


}; // end namespace

#endif /* TCALLER_H_ */
