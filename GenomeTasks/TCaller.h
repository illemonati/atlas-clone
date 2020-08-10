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
	std::string _callerName;
	std::string _filenameExtention;

	//lookup stuff
	TGenotypeMap _genoMap;
	BAM::TQualityMap _qualMap;
	VCF::TVCFInfoFields _VCFInfoFields;
	VCF::TVCFGenotypeFields _VCFGenotypeFields;
	TRandomGenerator* _randomGenerator;

	//output choices
	bool _printSitesWithNoData;
	bool _printAltIfHomoRef;
	bool _allowTriallelicSites;
	std::string _missingGenotype;

	//output file
	std::string _filename;
	gz::ogzstream _vcf;
	bool _vcfOpen;
	std::string _genotypeFormatString;

	//temp variables for calling
	std::string _calledGenotype;
	std::vector<Base> _genotypesWithHighestMetric;
	Base referenceBase;
	std::vector<Base> _altAlleles; //order of Base enums: A, C, G, T, N
	TBaseCounts _alleleCounts;
	bool _allelesCounted;

	//genotype prior
	bool _usesPrior;
	TGenotypeData* _genotypePrior; //for callers using a prior. Note: all callers accept priors, but may not use them.
	bool _priorSet;

	//functions regarding VCF file
	void _setAcceptableFields(VCF::TVCFFieldVector* fields, std::string tags);
	void _printField(VCF::TVCFFieldVector* fields, std::string tag);
	void _writeVCFHeader(const std::string & sampleName);

	//function to write info fields
	std::vector< std::string (TCaller::*)(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) > _VCFInfoFunctionsVec;
	void _fillInfoFieldFunctionPointers();
	virtual std::string _getVCFInfoString_DP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

	//functions to write genotype fields
	std::vector< std::string (TCaller::*)(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) > _VCFGenotypeFunctionsVec;
	void _fillGenotypeFieldFunctionPointers();
	virtual std::string _getVCFGenotypeString_GT(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	virtual std::string _getVCFGenotypeString_DP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	virtual std::string _getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw std::runtime_error("Function std::string getVCFGenotypeString_GQ(const const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };
	virtual std::string _getVCFGenotypeString_AD(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	virtual std::string _getVCFGenotypeString_AP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw std::runtime_error("Function std::string getVCFGenotypeString_AP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };
	virtual std::string _getVCFGenotypeString_GL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw std::runtime_error("Function std::string getVCFGenotypeString_GL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };
	virtual std::string _getVCFGenotypeString_PL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw std::runtime_error("Function std::string getVCFGenotypeString_PL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };
	virtual std::string _getVCFGenotypeString_GP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw std::runtime_error("Function std::string getVCFGenotypeString_GP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };
	virtual std::string _getVCFGenotypeString_AB(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw std::runtime_error("Function std::string getVCFGenotypeString_AB(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };
	virtual std::string _getVCFGenotypeString_AI(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){ throw std::runtime_error("Function std::string getVCFGenotypeString_AI(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };

	//write VCF
	std::string _composeVCFString(std::vector<std::string (TCaller::*)(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods)> & vec, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	virtual void _writeAlternativeAllelesToVCF();
	void _writeCallToVCF(const std::string & chr, const long pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void _writeMissingDataToVCF(const TSite & site);
	virtual void _clearAfterCall();

	//call
	void _countAlleles(const TSite & site);
	virtual void _callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	virtual void _callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

	template <typename T> uint8_t _pickIndexWithHighestMetric(T* metric, const uint8_t size){
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
		return vec[_randomGenerator->pickOne(vec.size())];
	};

	template <typename T> uint8_t _pickIndexWithSecondHighestMetric(T* metric, const uint8_t size, const uint8_t excludeIndex){
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
		return vec[_randomGenerator->pickOne(vec.size())];
	};

public:
	TCaller(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	virtual ~TCaller();

	std::string name() const{ return _callerName; };

	//set which fields to print
	void printInfoFields(std::vector<std::string> & tags);
	void printInfoFields(std::string tags);
	void printGenotypeFields(std::vector<std::string> & tags);
	void printGenotypeFields(std::string tags);

	//open / close _vcf file
	void openVCF(const std::string Filename, const std::string sampleName, TLog* logfile);
	void closeVCF();

	//other output options
	void initializeOutput(TParameters & Parameters, TLog* Logfile);
	bool printSitesWithNoData(){ return _printSitesWithNoData; };

	//prior
	bool usesPrior(){ return _usesPrior; };
	void setPrior(TGenotypeData* prior){ _genotypePrior = prior; _priorSet = true; };
	void call(const std::string & chr, const long pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void call(const std::string & chr, const long pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods, const char & firstAllele, const char & secondAllele);
};

//------------------------------------------------------
// TCallerRandomBase
//------------------------------------------------------
class TCallerRandomBase:public TCaller{
private:
	void _callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void _callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCallerRandomBase(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerMajorityCall
//------------------------------------------------------
class TCallerMajorityBase:public TCaller{
private:
	uint32_t _downsampleDepth;

	void _callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void _callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCallerMajorityBase(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
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
	void _callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void _callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string _getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string _getVCFGenotypeString_AP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCallerAllelePresence(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
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

	void _clearAfterCall();
	void callGenotypeFromMetric(TGenotypeData & metric);
	void callGenotypeFromMetricKnownAlleles(const TGenotypeData & metric, std::vector<int> & indeces);
	void callGenotypeFromMetricKnownAllelesUpdateIndex(const TGenotypeData & metric);
	std::string getPerGenotypeMetricString(TGenotypeData & metric);
	void calculateImbalance(const TSite & site);
	std::string _getVCFGenotypeString_AB(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string _getVCFGenotypeString_AI(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCallerDiploid(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerMLE
//------------------------------------------------------
class TCallerMLE:public TCallerDiploid{
private:
	void _callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void _callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string _getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string _getVCFGenotypeString_GL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string _getVCFGenotypeString_PL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCallerMLE(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerBayes
//------------------------------------------------------
class TCallerBayes:public TCallerDiploid{
private:
	TGenotypeProbabilities posterior;

	void _callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void _callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string _getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string _getVCFGenotypeString_GP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	std::string getVCFGenotypeString_PP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCallerBayes(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
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
