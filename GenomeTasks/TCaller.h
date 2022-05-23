/*
 * TCaller.h
 *
 *  Created on: Nov 17, 2018
 *      Author: phaentu
 */

#ifndef TCALLER_H_
#define TCALLER_H_

#include <stdexcept>
#include <stdint.h>
#include <string>
#include <vector>

#include "GenotypeTypes.h"
#include "TGenome.h"
#include "TGenotypeData.h"
#include "TTask.h"
#include "VCF/TVCFFields.h"
#include "gzstream.h"
#include "mathFunctions.h"
#include "stringFunctions.h"

namespace GenotypeLikelihoods { class TGenotypePrior; }
namespace GenotypeLikelihoods { class TSite; }

namespace GenomeTasks{



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
	genometools::TVCFInfoFields _VCFInfoFields;
	genometools::TVCFGenotypeFields _VCFGenotypeFields;

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
	GenotypeLikelihoods::TBaseCounts _alleleCounts;
	bool _allelesCounted;

	//genotype prior
	bool _usesPrior;
	GenotypeLikelihoods::TGenotypeProbabilities* _genotypePrior; //for callers using a prior. Note: all callers accept priors, but may not use them.
	bool _priorSet;

	//functions regarding VCF file
	void _setAcceptableFields(genometools::TVCFFieldVector* fields, std::string tags);
	void _printField(genometools::TVCFFieldVector* fields, std::string tag);
	void _writeVCFHeader(const std::string & sampleName);

	//function to write info fields
	std::vector< std::string (TCaller::*)(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) > _VCFInfoFunctionsVec;
	void _fillInfoFieldFunctionPointers();
	virtual std::string _getVCFInfoString_DP(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);

	//functions to write genotype fields
	std::vector< std::string (TCaller::*)(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) > _VCFGenotypeFunctionsVec;
	void _fillGenotypeFieldFunctionPointers();
	virtual std::string _getVCFGenotypeString_GT(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);
	virtual std::string _getVCFGenotypeString_DP(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);
	virtual std::string _getVCFGenotypeString_GQ(const GenotypeLikelihoods::TSite &, GenotypeLikelihoods::TGenotypeLikelihoods &){ throw std::runtime_error("Function std::string getVCFGenotypeString_GQ(const const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };
	virtual std::string _getVCFGenotypeString_AD(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);
	virtual std::string _getVCFGenotypeString_AP(const GenotypeLikelihoods::TSite &, GenotypeLikelihoods::TGenotypeLikelihoods &){ throw std::runtime_error("Function std::string getVCFGenotypeString_AP(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };
	virtual std::string _getVCFGenotypeString_GL(const GenotypeLikelihoods::TSite &, GenotypeLikelihoods::TGenotypeLikelihoods &){ throw std::runtime_error("Function std::string getVCFGenotypeString_GL(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };
	virtual std::string _getVCFGenotypeString_PL(const GenotypeLikelihoods::TSite &, GenotypeLikelihoods::TGenotypeLikelihoods &){ throw std::runtime_error("Function std::string getVCFGenotypeString_PL(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };
	virtual std::string _getVCFGenotypeString_GP(const GenotypeLikelihoods::TSite &, GenotypeLikelihoods::TGenotypeLikelihoods &){ throw std::runtime_error("Function std::string getVCFGenotypeString_GP(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };
	virtual std::string _getVCFGenotypeString_AB(const GenotypeLikelihoods::TSite &, GenotypeLikelihoods::TGenotypeLikelihoods &){ throw std::runtime_error("Function std::string getVCFGenotypeString_AB(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };
	virtual std::string _getVCFGenotypeString_AI(const GenotypeLikelihoods::TSite &, GenotypeLikelihoods::TGenotypeLikelihoods &){ throw std::runtime_error("Function std::string getVCFGenotypeString_AI(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) not defined for base class TCaller!"); };

	//write VCF
	std::string _composeVCFString(std::vector<std::string (TCaller::*)(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods)> & vec, const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);
	virtual void _writeAlternativeAllelesToVCF();
	void _writeCallToVCF(const std::string & chr, const long pos, const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);
	void _writeMissingDataToVCF(const GenotypeLikelihoods::TSite & site);
	virtual void _clearAfterCall();

	//call
	void _countAlleles(const GenotypeLikelihoods::TSite & site);
	virtual bool _callGenotype(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);
	virtual bool _callGenotypeKnownAlleles(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCaller();
	virtual ~TCaller();

	std::string name() const{ return _callerName; };

	//set which fields to print
	void printInfoFields(std::vector<std::string> & tags);
	void printInfoFields(std::string tags);
	void printGenotypeFields(std::vector<std::string> & tags);
	void printGenotypeFields(std::string tags);

	//open / close _vcf file
	void openVCF(const std::string Filename, const std::string sampleName);
	void closeVCF();

	//other output options
	void initializeOutput();
	bool printSitesWithNoData(){ return _printSitesWithNoData; };

	//prior
	bool usesPrior(){ return _usesPrior; };
	void setPrior(GenotypeLikelihoods::TGenotypeProbabilities* prior){ _genotypePrior = prior; _priorSet = true; };
	void call(const std::string & chr, long pos, const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);
	void call(const std::string & chr, long pos, const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods, const genometools::Base & firstAllele, const genometools::Base & secondAllele);
};

//------------------------------------------------------
// TCallerRandomBase
//------------------------------------------------------
class TCallerRandomBase:public TCaller{
private:
	bool _callGenotype(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	bool _callGenotypeKnownAlleles(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;

public:
	TCallerRandomBase();
};

//------------------------------------------------------
// TCallerMajorityCall
//------------------------------------------------------
class TCallerMajorityBase:public TCaller{
private:
	bool _callGenotype(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	bool _callGenotypeKnownAlleles(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;

public:
	TCallerMajorityBase();
};

//------------------------------------------------------
// TCallerConsensify
//------------------------------------------------------
class TCallerConsensify:public TCaller{
private:
	uint32_t _downsampleDepth, _minMajorityDepth;

	bool _callGenotype(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	void _callGenotypeKnownAlleles(const GenotypeLikelihoods::TBaseCounts & AlleleCounts);
	bool _callGenotypeKnownAlleles(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;

public:
	TCallerConsensify(uint32_t DownsampleDepth);
};

//------------------------------------------------------
// TCallerAllelePresence
//------------------------------------------------------
class TCallerAllelePresence:public TCaller{
private:
	GenotypeLikelihoods::TGenotypeProbabilities _posterior;
	GenotypeLikelihoods::TBaseLikelihoods _allelePostProb;
	genometools::Base _MAP;

	void _fillPosteriors(GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);
	bool _callGenotype(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	bool _callGenotypeKnownAlleles(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_GQ(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_AP(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;

public:
	TCallerAllelePresence();
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

	void _clearAfterCall() override;
	void callGenotypeFromMetric(const GenotypeLikelihoods::TGenotypeProbabilities & metric);
	void callGenotypeFromMetricKnownAlleles(const GenotypeLikelihoods::TGenotypeProbabilities & metric);
	bool callGenotypeFromMetricKnownAllelesUpdateIndex(const GenotypeLikelihoods::TGenotypeProbabilities & metric);

	template<typename T> std::string _getPerGenotypeMetricString(const T & metric){
		using genometools::Base;
		using coretools::str::toString;
		//if you have alleles R, A, B, C then the order is: RR, RA, AA | RB, AB, BB | RC, AC, BC, CC
		//plot missing value (.) for all metrics involving the reference if the reference is N
		std::string ret;
		//first for reference base
		if(referenceBase == Base::N)
			ret = ".";
		else
			ret = toString(metric[genometools::genotype(referenceBase, referenceBase)]);

		//now for alternative alleles
		if(_altAlleles.size() > 0){
			if(referenceBase == Base::N)
				ret += ",.";
			else
				ret += ',' + toString(metric[genometools::genotype(referenceBase, _altAlleles[0])]);
			ret += ',' + toString(metric[genometools::genotype(_altAlleles[0], _altAlleles[0])]);

			if(_altAlleles.size() > 1){
				if(referenceBase == Base::N)
					ret += ",.";
				else
					ret += ',' + toString(metric[genometools::genotype(referenceBase, _altAlleles[1])]);
				ret += ',' + toString(metric[genometools::genotype(_altAlleles[0], _altAlleles[1])]);
				ret += ',' + toString(metric[genometools::genotype(_altAlleles[1], _altAlleles[1])]);
			}

			if(_altAlleles.size() > 2){
				if(referenceBase == Base::N)
					ret += ",.";
				else
					ret += ',' + toString(metric[genometools::genotype(referenceBase, _altAlleles[2])]);
				ret += ',' + toString(metric[genometools::genotype(_altAlleles[0], _altAlleles[2])]);
				ret += ',' + toString(metric[genometools::genotype(_altAlleles[1], _altAlleles[2])]);
				ret += ',' + toString(metric[genometools::genotype(_altAlleles[2], _altAlleles[2])]);
			}
		}
		return ret;
	};

	void calculateImbalance(const GenotypeLikelihoods::TSite & site);
	std::string _getVCFGenotypeString_AB(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_AI(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;

public:
	TCallerDiploid();
};

//------------------------------------------------------
// TCallerMLE
//------------------------------------------------------
class TCallerMLE:public TCallerDiploid{
private:
	bool _callGenotype(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	bool _callGenotypeKnownAlleles(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_GQ(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_GL(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_PL(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;

public:
	TCallerMLE();
};

//------------------------------------------------------
// TCallerBayes
//------------------------------------------------------
class TCallerBayes:public TCallerDiploid{
private:
	GenotypeLikelihoods::TGenotypeProbabilities _posterior;

	bool _callGenotype(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	bool _callGenotypeKnownAlleles(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_GQ(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_GP(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods) override;
	std::string _getVCFGenotypeString_PP(const GenotypeLikelihoods::TSite & site, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);

public:
	TCallerBayes();
};

//------------------------------------------------------
// TCall
// the class to perform calls based on windows
//------------------------------------------------------
class TCall:public TGenome_windows{
private:
	TCaller* _caller;
	GenotypeLikelihoods::TGenotypePrior* _prior;

	void _initializeGenotypePrior();
	void _call();
	void _callKnwonAlleles();
	void _handleWindow() override;
	void _handleAlignment() override {}

public:
	TCall();
	~TCall();
	void call();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_call:public coretools::TTask{
public:
	TTask_call(){ _explanation = "Calling genotypes"; };

	void run(){
		using namespace coretools::instances;
		TCall caller;
		caller.call();
	};
};


}; // end namespace

#endif /* TCALLER_H_ */
