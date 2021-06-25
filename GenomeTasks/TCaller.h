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
	VCF::TVCFInfoFields _VCFInfoFields;
	VCF::TVCFGenotypeFields _VCFGenotypeFields;
	TRandomGenerator* _randomGenerator;

	//output choices
	bool _printSitesWithNoData;
	bool _printAltIfHomoRef;
	bool _allowTriallelicSites;
	bool _allowKnownAllelesCallsDifferentFromBestCall;
	std::string _missingGenotype;

	//output file
	std::string _filename;
	gz::ogzstream _vcf;
	bool _vcfOpen;
	std::string _genotypeFormatString;

	//temp variables for calling
	std::string _calledGenotype;
	std::vector<genometools::Base> _genotypesWithHighestMetric;
	genometools::Base referenceBase;
	std::vector<genometools::Base> _altAlleles;
	TBaseCounts _alleleCounts;
	bool _allelesCounted;

	//genotype prior
	bool _usesPrior;
	TGenotypeProbabilities* _genotypePrior; //for callers using a prior. Note: all callers accept priors, but may not use them.
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
	virtual bool _callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	virtual bool _callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

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
	void setPrior(TGenotypeProbabilities* prior){ _genotypePrior = prior; _priorSet = true; };
	void call(const std::string & chr, const long & pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);
	void call(const std::string & chr, const long & pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods, const genometools::Base & firstAllele, const genometools::Base & secondAllele);
};

//------------------------------------------------------
// TCallerRandomBase
//------------------------------------------------------
class TCallerRandomBase:public TCaller{
private:
	bool _callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	bool _callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;

public:
	TCallerRandomBase(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerMajorityCall
//------------------------------------------------------
class TCallerMajorityBase:public TCaller{
private:
	bool _callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	bool _callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;

public:
	TCallerMajorityBase(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerConsensify
//------------------------------------------------------
class TCallerConsensify:public TCaller{
private:
	uint32_t _downsampleDepth, _minMajorityDepth;

	bool _callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	void _callGenotypeKnownAlleles(const TBaseCounts & AlleleCounts);
	bool _callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;

public:
	TCallerConsensify(const uint32_t & DownsampleDepth, TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerAllelePresence
//------------------------------------------------------
class TCallerAllelePresence:public TCaller{
private:
	TGenotypeProbabilities posterior;
	TBaseLikelihoods allelePostProb;
	genometools::Base MAP;

	void _fillPosteriors(TGenotypeLikelihoods & genotypeLikelihoods);
	bool _callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	bool _callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_AP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;

public:
	TCallerAllelePresence(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerDiploid
//------------------------------------------------------
class TCallerDiploid:public TCaller{
protected:
	//uint8_t indexOfMax, indexOfSecond;
	genometools::Genotype genotypeAtMax, genotypeAtSecond;
	std::string AB, AI;
	bool imbalanceCalculated;
	//TGenotypeData tmpGenoData;
	coretools::TBinomPValue _binomP;

	void _clearAfterCall();
	void callGenotypeFromMetric(const TGenotypeProbability_base & metric);
	void callGenotypeFromMetricKnownAlleles(const TGenotypeProbability_base & metric);
	bool callGenotypeFromMetricKnownAllelesUpdateIndex(const TGenotypeProbability_base & metric);

	template<typename T> std::string _getPerGenotypeMetricString(const T & metric){
		//if you have alleles R, A, B, C then the order is: RR, RA, AA | RB, AB, BB | RC, AC, BC, CC
		//plot missing value (.) for all metrics involving the reference if the reference is N
		std::string ret;
		//first for reference base
		if(referenceBase == genometools::N)
			ret = ".";
		else
			ret = toString(metric[genometools::Genotype(referenceBase, referenceBase)]);

		//now for alternative alleles
		if(_altAlleles.size() > 0){
			if(referenceBase == genometools::N)
				ret += ",.";
			else
				ret += ',' + toString(metric[genometools::Genotype(referenceBase, _altAlleles[0])]);
			ret += ',' + toString(metric[genometools::Genotype(_altAlleles[0], _altAlleles[0])]);

			if(_altAlleles.size() > 1){
				if(referenceBase == genometools::N)
					ret += ",.";
				else
					ret += ',' + toString(metric[genometools::Genotype(referenceBase, _altAlleles[1])]);
				ret += ',' + toString(metric[genometools::Genotype(_altAlleles[0], _altAlleles[1])]);
				ret += ',' + toString(metric[genometools::Genotype(_altAlleles[1], _altAlleles[1])]);
			}

			if(_altAlleles.size() > 2){
				if(referenceBase == genometools::N)
					ret += ",.";
				else
					ret += ',' + toString(metric[genometools::Genotype(referenceBase, _altAlleles[2])]);
				ret += ',' + toString(metric[genometools::Genotype(_altAlleles[0], _altAlleles[2])]);
				ret += ',' + toString(metric[genometools::Genotype(_altAlleles[1], _altAlleles[2])]);
				ret += ',' + toString(metric[genometools::Genotype(_altAlleles[2], _altAlleles[2])]);
			}
		}
		return ret;
	};

	void calculateImbalance(const TSite & site);
	std::string _getVCFGenotypeString_AB(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_AI(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;

public:
	TCallerDiploid(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerMLE
//------------------------------------------------------
class TCallerMLE:public TCallerDiploid{
private:
	bool _callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	bool _callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_GL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_PL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;

public:
	TCallerMLE(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
};

//------------------------------------------------------
// TCallerBayes
//------------------------------------------------------
class TCallerBayes:public TCallerDiploid{
private:
	TGenotypeProbabilities posterior;

	bool _callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	bool _callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_GP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_PP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods);

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
class TTask_call:public coretools::TTask{
public:
	TTask_call(){ _explanation = "Calling genotypes"; };

	void run(coretools::TParameters & Parameters, coretools::TLog* Logfile){
		TCall caller(Parameters, Logfile, _randomGenerator);
		caller.call();
	};
};


}; // end namespace

#endif /* TCALLER_H_ */
