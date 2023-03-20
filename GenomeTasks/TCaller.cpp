/*
 * TCaller.cpp
 *
 *  Created on: Nov 17, 2018
 *      Author: phaentu
 */

#include "TCaller.h"

#include <math.h>
#include <stddef.h>
#include <algorithm>
#include <exception>
#include <memory>
#include <ostream>
#include <set>

#include "coretools/Main/TError.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "TGenotypeData.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TGenotypePrior.h"
#include "coretools/Containers/TMassFunction.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TSequencedBase.h"
#include "TSite.h"
#include "TSiteSubset.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Containers/TMassFunction.h"
#include "genometools/VCF/TVCFFields.h"
#include "TWindow.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Types/probability.h"
#include "coretools/Types/weakTypes.h"

namespace GenomeTasks{

using namespace GenotypeLikelihoods;
using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using coretools::Probability;
using genometools::Base;
using namespace coretools::str;

namespace /*anonymous*/ {
template<template<typename, typename, size_t, typename...> typename Container, typename Type, typename Index, size_t N, typename... Args>
auto sampleFirstSecond(const Container<Type, Index, N, Args...> &c) {
	std::array<Index, N> is;
	for (size_t i = 0; i < N; ++i) is[i] = Index(i);
	std::sort(is.begin(), is.end(), [&c](auto i, auto j){return c[i] > c[j];});

	size_t count_1  = 0;
	const auto max1 = c[is.front()];
	while (c[is[count_1]] == max1) {++count_1;}

	if(count_1 == 1){
		size_t count_2  = 1;
		const auto max2 = c[is[1]];
		while (c[is[count_2]] == max2) {++count_2;}

		return std::make_tuple(is[0], is[randomGenerator().sample(count_2) + 1]);
	} else if(count_1 == 2){
		const auto i1 = randomGenerator().sample(count_1);
		return std::make_tuple(is[i1], is[1-i1]);
	} else { // more than two max
		const auto i1 = randomGenerator().sample(count_1);
		auto i2 = i1;
		while (i2 == i1) i2 = randomGenerator().sample(count_1);
		return std::make_tuple(is[i1], is[i2]);
	}
}
} // namespace

/////////////////////////////////////////////////////////
// TCaller
/////////////////////////////////////////////////////////
TCaller::TCaller(){
	_callerName = "default caller";
	_filenameExtention = ".vcf.gz";

	//output choices
	_printSitesWithNoData = false;
	_printAltIfHomoRef = true;
	_allowTriallelicSites = true;
	_allowKnownAllelesCallsDifferentFromBestCall = true;
	_missingGenotype = ".";

	//vcf file
	_vcfOpen = false;

	//set acceptable tags
	_setAcceptableFields(&_VCFInfoFields, "DP");
	_setAcceptableFields(&_VCFGenotypeFields, "GT,DP,AD");

	//prior
	_genotypePrior = nullptr;
	_usesPrior = false;
	_priorSet = false;

	//set default tags to print
	printInfoFields("");
	printGenotypeFields("GT,DP");

	//tmp variables
	_allelesCounted = false;
};

TCaller::~TCaller(){
	closeVCF();
};

//-------------------------------------------------------------------------------------------
// Output settings
//-------------------------------------------------------------------------------------------
void TCaller::_setAcceptableFields(genometools::TVCFFieldVector* fields, std::string tags){
	std::vector<std::string> vec;
	fillContainerFromStringAny(tags, vec, ",", true);
	for(std::vector<std::string>::iterator it = vec.begin(); it != vec.end(); ++it)
		fields->acceptField(*it);
};

void TCaller::_printField(genometools::TVCFFieldVector* fields, std::string tag){
	if(!fields->useField(tag))
		UERROR("VCF ", fields->type(), " field '", tag, "' can not be printed by the ", _callerName, "!");
};

void TCaller::printInfoFields(std::vector<std::string> & tags){
	_VCFInfoFields.clearUsed();
	for(std::vector<std::string>::iterator it = tags.begin(); it != tags.end(); ++it){
		_printField(&_VCFInfoFields, *it);
	}
	_fillInfoFieldFunctionPointers();
};

void TCaller::printInfoFields(std::string tags){
	std::vector<std::string> vec;
	fillContainerFromStringAny(tags, vec, ",", true);
	printInfoFields(vec);
};

void TCaller::printGenotypeFields(std::vector<std::string> & tags){
	_VCFGenotypeFields.clearUsed();
	for(std::vector<std::string>::iterator it = tags.begin(); it != tags.end(); ++it)
		_printField(&_VCFGenotypeFields, *it);
	_fillGenotypeFieldFunctionPointers();
};

void TCaller::printGenotypeFields(std::string tags){
	std::vector<std::string> vec;
	fillContainerFromStringAny(tags, vec, ",", true);
	printGenotypeFields(vec);
};

void TCaller::initializeOutput(){
	//info fields
	if(parameters().parameterExists("infoFields")){
		logfile().listFlush("Parsing VCF info fields ...");
		printInfoFields(parameters().getParameter<std::string>("infoFields"));
		logfile().done();
	}
	if(_VCFInfoFields.numUsed() > 0){
		logfile().list("Will print these VCF info fields: " + _VCFInfoFields.getListOfUsedFields(", ") + ". (parameter 'infoFields')");
	} else {
		logfile().list("Will not print any VCF info fields. (parameter 'infoFields')");
	}

	//genotype fields
	if(parameters().parameterExists("formatFields")){
		logfile().listFlush("Parsing VCF format fields ...");
		printGenotypeFields(parameters().getParameter<std::string>("formatFields"));
		logfile().done();
	}
	if(_VCFGenotypeFields.numUsed() > 0){
		logfile().list("Will print these VCF format fields: " + _VCFGenotypeFields.getListOfUsedFields(", ") + ". (parameter 'formatFields')");
	} else {
		logfile().list("Will not print any VCF format fields. (parameter 'formatFields')");
	}

	//other output options
	if(parameters().parameterExists("printAll")){
		_printSitesWithNoData = true;
		logfile().list("Will print all sites, also those without data. (parameter 'printAll')");
	} else {
		_printSitesWithNoData = false;
		logfile().list("Will print only sites with data. (use 'printAll' to print all);");
	}

	if(parameters().parameterExists("noAltIfHomoRef")){
		_printAltIfHomoRef = false;
		logfile().list("Will not print an alternative allele if the call is homozygous reference. (parameter 'noAltIfHomoRef')");
	} else {
		_printAltIfHomoRef = true;
		logfile().list("Will print the most likely alternative allele even if the call is homozygous reference. (use 'noAltIfHomoRef' to suppress)");
	}

	if(parameters().parameterExists("noTriallelic")){
		_allowTriallelicSites = false;
		logfile().list("Will not call genotypes resulting in two alternative alleles. (parameter 'noTriallelic')");
	} else {
		_allowTriallelicSites = true;
		logfile().list("Will allow for genotypes with two alternative alleles. (use 'noTriallelic' to suppress)");
	}

	if(parameters().parameterExists("noCallsViolatingBest")){
		_allowKnownAllelesCallsDifferentFromBestCall = false;
		logfile().list("Will not call genotypes from known alleles that conflict with best call across all genotypes. (parameter 'noCallsViolatingBest')");
	} else {
		_allowKnownAllelesCallsDifferentFromBestCall = true;
		logfile().list("Will call genotypes from known alleles even if they differ from best call across all genotypes. (use 'noCallsViolatingBest' to suppress)");
	}
}

//-------------------------------------------------------------------------------------------
// open / close VCF file, print header
//-------------------------------------------------------------------------------------------
void TCaller::openVCF(const std::string FilenameTag, const std::string sampleName){
	_filename = FilenameTag  + _filenameExtention + ".gz";
	logfile().list("Writing calls to VCF file '" + _filename + "'.");
	_vcf.open(_filename.c_str());
	if(!_vcf) UERROR("Failed to open VCF file '", _filename, "' for writing!");
	_vcfOpen = true;

	//write header
	_writeVCFHeader(sampleName);
};

void TCaller::closeVCF(){
	if(_vcfOpen){
		_vcf.close();
		_vcfOpen = false;
	}
};

void TCaller::_writeVCFHeader(const std::string & sampleName){
	//write header
	_vcf << "##fileformat=VCFv4.2\n";
	_vcf << "##source=atlas\n";

	//write INFO and GENOTYPE fields
	_VCFInfoFields.writeVCFHeader(_vcf);
	_VCFGenotypeFields.writeVCFHeader(_vcf);

	//write column header
	_vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleName << "\n";
};

//-------------------------------------------------------------------------------------------
// Info fields
//-------------------------------------------------------------------------------------------
void TCaller::_fillInfoFieldFunctionPointers(){
	//clear current vector
	_VCFInfoFunctionsVec.clear();

	//get used tags
	std::vector<std::string> tagVec;
	_VCFInfoFields.fillVectorWithTagsOfUsedFields(tagVec);

	//add functions to info field vector
	for(std::vector<std::string>::iterator it = tagVec.begin(); it != tagVec.end(); ++it){
		if(*it == "DP")
			_VCFInfoFunctionsVec.push_back( &TCaller::_getVCFInfoString_DP );
		else UERROR("No function defined for VCF ", _VCFInfoFields.type(), " field '", *it, "'! @Programmer: add function to TTCaller::fillInfoFieldFunctionPointers()!");
	}

};

std::string TCaller::_getVCFInfoString_DP(const TSite & site, TGenotypeLikelihoods &){
	return "DP=" + toString(site.depth());
};

//-------------------------------------------------------------------------------------------
// genotype fields
//-------------------------------------------------------------------------------------------
void TCaller::_fillGenotypeFieldFunctionPointers(){
	//clear current vector
	_VCFGenotypeFunctionsVec.clear();
	_genotypeFormatString.clear();

	//get used tags
	std::vector<std::string> tagVec;
	_VCFGenotypeFields.fillVectorWithTagsOfUsedFields(tagVec);

	//add functions to genotype field vector
	for(std::vector<std::string>::iterator it = tagVec.begin(); it != tagVec.end(); ++it){
		if(*it == "GT")
			_VCFGenotypeFunctionsVec.push_back( &TCaller::_getVCFGenotypeString_GT );
		else if(*it == "DP")
			_VCFGenotypeFunctionsVec.push_back( &TCaller::_getVCFGenotypeString_DP );
		else if(*it == "GQ")
			_VCFGenotypeFunctionsVec.push_back( &TCaller::_getVCFGenotypeString_GQ );
		else if(*it == "AD")
			_VCFGenotypeFunctionsVec.push_back( &TCaller::_getVCFGenotypeString_AD );
		else if(*it == "AP")
			_VCFGenotypeFunctionsVec.push_back( &TCaller::_getVCFGenotypeString_AP );
		else if(*it == "GL")
			_VCFGenotypeFunctionsVec.push_back( &TCaller::_getVCFGenotypeString_GL );
		else if(*it == "PL")
			_VCFGenotypeFunctionsVec.push_back( &TCaller::_getVCFGenotypeString_PL );
		else if(*it == "GP")
			_VCFGenotypeFunctionsVec.push_back( &TCaller::_getVCFGenotypeString_GP );
		else if(*it == "AB")
			_VCFGenotypeFunctionsVec.push_back( &TCaller::_getVCFGenotypeString_AB );
		else if(*it == "AI")
			_VCFGenotypeFunctionsVec.push_back( &TCaller::_getVCFGenotypeString_AI );
		else UERROR("No function defined for VCF ", _VCFGenotypeFields.type(), " field '", *it, "'! @Programmer: add function to TTCaller::fillGenotypeFieldFunctionPointers()!");

		//add to format string
		if(_genotypeFormatString.length() > 0) _genotypeFormatString += ':';
		_genotypeFormatString += *it;
	}
};

std::string TCaller::_getVCFGenotypeString_GT(const TSite &, TGenotypeLikelihoods &){
	return _calledGenotype;
};

std::string TCaller::_getVCFGenotypeString_DP(const TSite & site, TGenotypeLikelihoods &){
	return toString(site.depth());
};

std::string TCaller::_getVCFGenotypeString_AD(const TSite & site, TGenotypeLikelihoods &){
	_countAlleles(site);
	std::string ret;
	if(referenceBase == genometools::Base::N) ret = "0";
	else ret = toString(_alleleCounts[referenceBase]);

	for(auto& a : _altAlleles)
		ret += ',' + toString(_alleleCounts[a]);
	return ret;
};

//-------------------------------------------------------------------------------------------
// writing VCF
//-------------------------------------------------------------------------------------------
std::string TCaller::_composeVCFString(std::vector<std::string (TCaller::*)(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods)> & vec, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//no info fields?
	if(vec.empty()) return ".";

	//add first info
	auto it = vec.begin();
	std::string info = (this->*(*it))(site, genotypeLikelihoods);
	++it;

	//loop over rest
	for(; it != vec.end(); ++it)
		info += ':' + (this->*(*it))(site, genotypeLikelihoods);

	return info;
};

void TCaller::_writeAlternativeAllelesToVCF(){
	if(_altAlleles.size() == 0){
		_vcf << '.';
	} else {
		_vcf << _altAlleles[0];
		for(size_t i=1; i<_altAlleles.size(); ++i)
			_vcf << ',' << _altAlleles[i];
	}
};

void TCaller::_writeCallToVCF(const std::string & chr, const long pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//apply filter on alternative alleles
	if(!_printAltIfHomoRef && (_calledGenotype == "0/0" || _calledGenotype == "0"))
		_altAlleles.clear();

	//write chr, position and (no) variant ID
	_vcf << chr << '\t' << pos + 1 << "\t.\t"; //all internal positions are zero-based!

	//write reference and alternative alleles
	_vcf << site.refBase << "\t";
	_writeAlternativeAllelesToVCF();

	//write (no) variant quality and (no) filter
	_vcf << "\t.\t.";

	//write info fields
	_vcf << '\t' << _composeVCFString(_VCFInfoFunctionsVec, site, genotypeLikelihoods);

	//write genotype fields
	_vcf << '\t' << _genotypeFormatString << '\t' << _composeVCFString(_VCFGenotypeFunctionsVec, site, genotypeLikelihoods);

	//end with new line
	_vcf << '\n';

	//clean up storage
	_clearAfterCall();
};

void TCaller::_writeMissingDataToVCF(const TSite & site){
	if(_printSitesWithNoData)
		_vcf << "\t.\t" << site.refBase << "\t.\t.\t.\t.\tGT:DP\t" << _missingGenotype << ":0";
};

void TCaller::_clearAfterCall(){
	_altAlleles.clear();
	_allelesCounted = false;
};

//-------------------------------------------------------------------------------------------
// calling
//-------------------------------------------------------------------------------------------
void TCaller::_countAlleles(const TSite & site){
	if(!_allelesCounted){
		site.countAlleles(_alleleCounts);
		_allelesCounted = true;
	}
};

bool TCaller::_callGenotype(const TSite &, TGenotypeLikelihoods &){
	_calledGenotype = "./.";
	return true;
};

bool TCaller::_callGenotypeKnownAlleles(const TSite &, TGenotypeLikelihoods &){
	_calledGenotype = "./.";
	return true;
};

void TCaller::call(const std::string & chr, long pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//set reference base from site
	referenceBase = site.refBase;

	//check if there is data
	if(site.empty() || (referenceBase == genometools::Base::N && !_allowTriallelicSites) || !_callGenotype(site, genotypeLikelihoods))
		_writeMissingDataToVCF(site);
	else {
		_writeCallToVCF(chr, pos, site, genotypeLikelihoods);
	}
};

void TCaller::call(const std::string & chr, long pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods, const genometools::Base & firstAllele, const genometools::Base & secondAllele){
	//check if there is data
	if(!site.empty()){
		//set reference base from site
		referenceBase = site.refBase;

		//call
		if(referenceBase == firstAllele)
			_altAlleles.push_back(secondAllele);
		else
			_altAlleles.push_back(firstAllele);

		if(!_callGenotypeKnownAlleles(site, genotypeLikelihoods)){
			_writeMissingDataToVCF(site);
		} else {
			_writeCallToVCF(chr, pos, site, genotypeLikelihoods);
		}

	} else {
		_writeMissingDataToVCF(site);
	}
};

/////////////////////////////////////////////////////////
// TCallerRandomBase
/////////////////////////////////////////////////////////
TCallerRandomBase::TCallerRandomBase(){
	//caller settings
	_callerName = "Random Base Caller";
	_filenameExtention = "_randomBase.vcf";
	logfile().list("Will use the " + _callerName + ".");

	//parse VCF fields
	_setAcceptableFields(&_VCFInfoFields, "DP");
	_setAcceptableFields(&_VCFGenotypeFields, "GT,DP,AD");
	initializeOutput();

};

bool TCallerRandomBase::_callGenotype(const TSite & site, TGenotypeLikelihoods &){
	//randomly pick a base
	genometools::Base allele = site[randomGenerator().sample(site.depth())].base;

	//decide on alt
	if(allele == referenceBase){
		_calledGenotype = "0";
	} else {
		_altAlleles.push_back(allele);
		_calledGenotype = "1";
	}

	return true;
};

bool TCallerRandomBase::_callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods &){
	if(_allowKnownAllelesCallsDifferentFromBestCall){
		//pick among known alleles only
		_countAlleles(site);
		double probRef = (double) _alleleCounts[referenceBase] / (double) (_alleleCounts[referenceBase] + _alleleCounts[_altAlleles[0]]);

		//pick among known alleles
		if(randomGenerator().getRand() < probRef){
			_calledGenotype = "0";
		} else {
			_calledGenotype = "1";
		}
		return true;
	} else {
		//pick among all alleles and check
		genometools::Base allele = site[randomGenerator().sample(site.depth())].base;
		if(allele == referenceBase){
			_calledGenotype = "0";
			return true;
		} else if(allele == _altAlleles[0]){
			_calledGenotype = "1";
			return true;
		} else {
			return false;
		}
	}
};

/////////////////////////////////////////////////////////
// TCallerMajorityCall
/////////////////////////////////////////////////////////
TCallerMajorityBase::TCallerMajorityBase(){
	//general caller settings
	_callerName = "Majority Base Caller";
	_filenameExtention = "_majorityBase.vcf";
	logfile().list("Will use the " + _callerName + ".");

	//parse VCF fields
	initializeOutput();
};

bool TCallerMajorityBase::_callGenotype(const TSite & site, TGenotypeLikelihoods &){
	//get per allele counts
	_countAlleles(site);

	//call majority
	const auto [majorityBase, second] = sampleFirstSecond(_alleleCounts);

	//decide on alt
	if(majorityBase == referenceBase){
		_calledGenotype = "0";

		//find second most common as alternative allele
		_altAlleles.push_back(second);
	} else {
		_altAlleles.push_back(majorityBase);
		_calledGenotype = "1";
	}

	return true;
};

bool TCallerMajorityBase::_callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods &){
	//get per allele counts
	_countAlleles(site);

	if(_allowKnownAllelesCallsDifferentFromBestCall){
		//pick among known alleles only
		if(_alleleCounts[referenceBase] > _alleleCounts[_altAlleles[0]]){
			_calledGenotype = "0";
		} else if(_alleleCounts[referenceBase] < _alleleCounts[_altAlleles[0]]){
			_calledGenotype = "1";
		} else {
			//equal counts: pick at random
			if(randomGenerator().getRand() < 0.5){
				_calledGenotype = "0";
			} else {
				_calledGenotype = "1";
			}
		}
		return true;
	} else {
		//pick among all alleles and check
		//call majority
		const auto majorityBase = randomGenerator().sampleIndexOfMaxima<TBaseCounts, Base, 4>(_alleleCounts);

		//decide on call
		if(majorityBase == referenceBase){
			_calledGenotype = "0";
			return true;
		} else if(majorityBase == _altAlleles[0]){
			_calledGenotype = "1";
			return true;
		}
		return false;
	}
};

/////////////////////////////////////////////////////////
// TCallerConsensify
/////////////////////////////////////////////////////////
TCallerConsensify::TCallerConsensify(uint32_t DownsampleDepth) : TCaller() {
	//general caller settings
	_callerName = "Consensify Base Caller";
	_filenameExtention = "_consensifyBase.vcf";
	logfile().list("Will use the " + _callerName + ".");

	//parse VCF fields
	initializeOutput();

	//check downsample depth
	if(DownsampleDepth < 1){
		UERROR("Consensify caller requires downsampling of reads! Use 'downsample' to specify.");
	}
	_downsampleDepth = DownsampleDepth;
	_minMajorityDepth = ceil(_downsampleDepth / 2.0);

	logfile().list("Will call based on a majority out of " + toString(_downsampleDepth) + " bases. (parameter 'downsample')");
};

bool TCallerConsensify::_callGenotype(const TSite & site, TGenotypeLikelihoods &){
	if(site.depth() > _downsampleDepth){
		DEVERROR("depth > _downsampleDepth!");
	}

	//get per allele counts
	_countAlleles(site);

	//call majority
	const auto [majorityBase, second] = sampleFirstSecond(_alleleCounts);

	//check if we have sufficient depth to call
	if(_alleleCounts[majorityBase] < _minMajorityDepth){
		return false;
	}

	//decide on alt
	if(majorityBase == referenceBase){
		_calledGenotype = "0";

		//find second most common as alternative allele
		_altAlleles.emplace_back(second);
	} else {
		_altAlleles.emplace_back(majorityBase);
		_calledGenotype = "1";
	}

	return true;
};

bool TCallerConsensify::_callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods &){
	if(site.depth() > _downsampleDepth){
		DEVERROR("depth > _downsampleDepth!");
	}

	//get per allele counts
	_countAlleles(site);

	// call majority
	const auto majorityBase = randomGenerator().sampleIndexOfMaxima<TBaseCounts, Base, 4>(_alleleCounts);

	//check if we have sufficient depth to call
	if(_alleleCounts[majorityBase] < _minMajorityDepth){
		return false;
	}

	//check if allele fits
	//NOTE: _allowKnownAllelesCallsDifferentFromBestCall has no effect
	if(majorityBase == referenceBase){
		_calledGenotype = "0";
		return true;
	} else if(majorityBase == _altAlleles[0]){
		_calledGenotype = "1";
		return true;
	} else {
		return false;
	}
};

/////////////////////////////////////////////////////////
// TCallerAllelePresence
/////////////////////////////////////////////////////////
TCallerAllelePresence::TCallerAllelePresence() : TCaller(){
	//caller settings
	_callerName = "Allele Presence Caller";
	_filenameExtention = "_allelePresence.vcf";
	logfile().list("Will use the " + _callerName + ".");
	_usesPrior = true;

	//parse VCF fields
	_setAcceptableFields(&_VCFGenotypeFields, "GT,DP,AD,GQ,AP");
	initializeOutput();

	//initialize allele counts
	_MAP = genometools::Base::N;
};

void TCallerAllelePresence::_fillPosteriors(TGenotypeLikelihoods & genotypeLikelihoods){
	using namespace genometools;
	//calculate posterior probabilities
	_posterior = posterior(genotypeLikelihoods, *_genotypePrior);

	//sum for each base
	_allelePostProb[Base::A] = _posterior[Genotype::AA] + _posterior[Genotype::AC] + _posterior[Genotype::AG] + _posterior[Genotype::AT];
	_allelePostProb[Base::C] = _posterior[Genotype::AC] + _posterior[Genotype::CC] + _posterior[Genotype::CG] + _posterior[Genotype::CT];
	_allelePostProb[Base::G] = _posterior[Genotype::AG] + _posterior[Genotype::CG] + _posterior[Genotype::GG] + _posterior[Genotype::GT];
	_allelePostProb[Base::T] = _posterior[Genotype::AT] + _posterior[Genotype::CT] + _posterior[Genotype::GT] + _posterior[Genotype::TT];
};

bool TCallerAllelePresence::_callGenotype(const TSite &, TGenotypeLikelihoods & genotypeLikelihoods){
	if(!_priorSet) UERROR("Can not call AllelePresence genotypes: prior has not been set!");

	//fill posteriors for each allele
	_fillPosteriors(genotypeLikelihoods);

	//find MAP
	const auto [majorityBase, second] = sampleFirstSecond(_alleleCounts);
	_MAP = majorityBase;

	//decide on alt
	if(_MAP == referenceBase){
		_calledGenotype = "0";

		//find second most common as alternative allele
		_altAlleles.push_back(second);
	} else {
		_altAlleles.push_back(_MAP);
		_calledGenotype = "1";
	}

	return true;
};

bool TCallerAllelePresence::_callGenotypeKnownAlleles(const TSite &, TGenotypeLikelihoods & genotypeLikelihoods){
	if(!_priorSet) UERROR("Can not call AllelePresence genotypes: prior has not been set!");

	//fill posteriors for each allele
	_fillPosteriors(genotypeLikelihoods);

	//find MAP
	if(_allelePostProb[_altAlleles[0]] > _allelePostProb[referenceBase]){
		_MAP = _altAlleles[0];
		_calledGenotype = "1";
	} else if(_allelePostProb[referenceBase] > _allelePostProb[_altAlleles[0]]){
		_MAP = referenceBase;
		_calledGenotype = "0";
	} else {
		//both are equal, pick at random
		if(randomGenerator().getRand() < 0.5){
			_MAP = _altAlleles[0];
			_calledGenotype = "1";
		} else {
			_MAP = referenceBase;
			_calledGenotype = "0";
		}
	}

	if(!_allowKnownAllelesCallsDifferentFromBestCall){
		if(_allelePostProb[_MAP] < *std::max_element(_allelePostProb.begin(), _allelePostProb.end())){
			return false;
		}
	}

	return true;
};

std::string TCallerAllelePresence::_getVCFGenotypeString_GQ(const TSite &, TGenotypeLikelihoods &){
	return (std::string) genometools::PhredIntProbability(_allelePostProb[_MAP].complement());
};

std::string TCallerAllelePresence::_getVCFGenotypeString_AP(const TSite &, TGenotypeLikelihoods &){
	using genometools::Base;
	std::string ret = (std::string) genometools::PhredIntProbability(_allelePostProb[Base::A]);
	ret += ',' + (std::string) genometools::PhredIntProbability(_allelePostProb[Base::C]);
	ret += ',' + (std::string) genometools::PhredIntProbability(_allelePostProb[Base::G]);
	ret += ',' + (std::string) genometools::PhredIntProbability(_allelePostProb[Base::T]);
	return ret;
};

/////////////////////////////////////////////////////////
// TCallerDiploid
// common function between MLE, Bayes and gvcf callers
/////////////////////////////////////////////////////////
TCallerDiploid::TCallerDiploid() : TCaller(){
	imbalanceCalculated = false;
	_missingGenotype = "./.";

	_setAcceptableFields(&_VCFGenotypeFields, "AB,AI");
};

void TCallerDiploid::_clearAfterCall(){
	TCaller::_clearAfterCall();
	imbalanceCalculated = false;
};

void TCallerDiploid::callGenotypeFromMetric(const TGenotypeProbabilities & metric){
	using namespace genometools;
	const auto [genotypeAtMax, genotypeAtSecond] = sampleFirstSecond(metric);

	//decide on alternative alleles
	if(first(genotypeAtMax) == referenceBase){
		if(second(genotypeAtMax) == referenceBase){
			_calledGenotype = "0/0";
			//MLE is homozygous reference -> find second best allele
			if(first(genotypeAtSecond) == referenceBase)
				_altAlleles.push_back(second(genotypeAtSecond));
			else if(second(genotypeAtSecond) == referenceBase)
				_altAlleles.push_back(first(genotypeAtSecond));
			else {
				//none is ref, pick at random
				const auto b = randomGenerator().getRand() < 0.5 ? first(genotypeAtSecond) : second(genotypeAtSecond);
				_altAlleles.push_back( b );
			}
		} else {
			_altAlleles.push_back(second(genotypeAtMax));
			_calledGenotype = "0/1A";
		}
	} else {
		if(second(genotypeAtMax) == referenceBase){
			_altAlleles.push_back(first(genotypeAtMax));
			_calledGenotype = "0/1B";
		} else {
			if(isHomozygous(genotypeAtMax)){
				_altAlleles.push_back(first(genotypeAtMax));
				_calledGenotype = "1/1";

				//find second best allele, but give preference to reference in case likelihoods are equal
				if(_allowTriallelicSites){
					//int hetRef = genoMap.getGenotype(referenceBase, genoMap.genotypeToBase[indexOfMax][0]);
					//int homRef = genoMap.getGenotype(referenceBase, referenceBase);

					//only use second alternative allele in case het genotype with reference is less likely
					if(referenceBase == Base::N ||
					   (metric[genotype(referenceBase, first(genotypeAtMax))] < metric[ genotypeAtSecond]
						&& metric[genotype(referenceBase, referenceBase) ] < metric[ genotypeAtSecond ])){
						if(first(genotypeAtSecond) == referenceBase || first(genotypeAtSecond) == _altAlleles[0])
							_altAlleles.push_back(second(genotypeAtSecond));
						else if(second(genotypeAtSecond) == referenceBase || second(genotypeAtSecond) == _altAlleles[0])
							_altAlleles.push_back(first(genotypeAtSecond));
						else {
							//none is ref, pick at random
							const auto b = randomGenerator().getRand() < 0.5 ? first(genotypeAtSecond) : second(genotypeAtSecond);
							_altAlleles.push_back(b);
						}
					}
				}
			} else {
				if(_allowTriallelicSites){
					//allow triallelic sites
					_altAlleles.push_back(first(genotypeAtMax));
					_altAlleles.push_back(second(genotypeAtMax));
					_calledGenotype = "1/2";
				} else {
					//decide on which of the two alternative alleles to pick -> check second highest
					if(isHomozygous(genotypeAtSecond)){
						if(first(genotypeAtSecond) == first(genotypeAtMax))
							_altAlleles.push_back(first(genotypeAtMax));
						else if(first(genotypeAtSecond) == second(genotypeAtMax))
							_altAlleles.push_back(second(genotypeAtMax));
						else {
							//neither alt 0 nor 1, pick at random
							const auto b = randomGenerator().getRand() < 0.5 ? first(genotypeAtMax) : second(genotypeAtMax);
							_altAlleles.push_back(b);
						}
					} else {
						if(first(genotypeAtSecond) == referenceBase && (second(genotypeAtSecond) == first(genotypeAtMax) || second(genotypeAtSecond) == first(genotypeAtMax))){
							_altAlleles.push_back(second(genotypeAtSecond));
						} else if(second(genotypeAtSecond) == referenceBase && (first(genotypeAtSecond) == first(genotypeAtMax) || first(genotypeAtSecond) == second(genotypeAtMax))){
							_altAlleles.push_back(first(genotypeAtSecond));
						} else {
							//neither alt 0 nor 1, pick at random
							const auto b = randomGenerator().getRand() < 0.5 ? first(genotypeAtMax) : second(genotypeAtMax);
							_altAlleles.push_back(b);
						}
					}

					//now call genotype from these alleles
					callGenotypeFromMetricKnownAlleles(metric);
				}
			}
		}
	}
};

void TCallerDiploid::callGenotypeFromMetricKnownAlleles(const TGenotypeProbabilities & metric){
	//get genotypes
	genometools::BiallelicGenotypesWithAlleles geno(referenceBase, _altAlleles[0]);

	//find max
	Probability max = metric[geno.homoFirst()];
	if(metric[geno.het()] > max) max = metric[geno.het()];
	if(metric[geno.homoSecond()] > max) max = metric[geno.homoSecond()];

	//fill vector of all
	std::vector<std::string> vec;
	if(metric[geno.homoFirst()] == max){
		vec.push_back("0/0");
	}
	if(metric[geno.het()] == max){
		vec.push_back("0/1");
	}
	if(metric[geno.homoSecond()] == max){
		vec.push_back("1/1");
	}

	_calledGenotype = vec[randomGenerator().sample(vec.size())];
};

bool TCallerDiploid::callGenotypeFromMetricKnownAllelesUpdateIndex(const TGenotypeProbabilities & metric){
	//get genotypes
	genometools::BiallelicGenotypesWithAlleles geno(referenceBase, _altAlleles[0]);

	//identify highest and second highest
	//TODO: is there a better way than an endless if / else?
	if(metric[geno.homoFirst()] > metric[geno.het()]){
		if(metric[geno.homoFirst()] > metric[geno.homoSecond()]){
			genotypeAtMax = geno.homoFirst();
			_calledGenotype = "0/0";
			if(metric[geno.het()] > metric[geno.homoSecond()]){
				genotypeAtSecond = geno.het();
			} else if(metric[geno.het()] < metric[geno.homoSecond()]){
				genotypeAtSecond = geno.homoSecond();
			} else {
				//het and homoSecond are equal: pick at random
				if(randomGenerator().getRand() < 0.5){
					genotypeAtSecond = geno.het();
				} else {
					genotypeAtSecond = geno.homoSecond();
				}
			}
		} else if(metric[geno.homoFirst()] < metric[geno.homoSecond()]){
			genotypeAtMax = geno.homoSecond();
			_calledGenotype = "1/1";
			genotypeAtSecond = geno.homoFirst();
		} else {
			//homoFirst and homoSecond are equal: pick at random
			if(randomGenerator().getRand() < 0.5){
				genotypeAtMax = geno.homoFirst();
				_calledGenotype = "0/0";
				genotypeAtSecond = geno.homoSecond();
			} else {
				genotypeAtMax = geno.homoSecond();
				_calledGenotype = "1/1";
				genotypeAtSecond = geno.homoFirst();
			}
		}
	} else if(metric[geno.homoFirst()] < metric[geno.het()]){
		if(metric[geno.het()] > metric[geno.homoSecond()]){
			genotypeAtMax = geno.het();
			_calledGenotype = "0/1";
			if(metric[geno.homoFirst()] > metric[geno.homoSecond()]){
				genotypeAtSecond = geno.homoFirst();
			} else if(metric[geno.homoFirst()] < metric[geno.homoSecond()]){
				genotypeAtSecond = geno.homoSecond();
			} else {
				//homoFirst and homoSecond are equal: pick at random
				if(randomGenerator().getRand() < 0.5){
					genotypeAtSecond = geno.homoFirst();
				} else {
					genotypeAtSecond = geno.homoSecond();
				}
			}
		} else if(metric[geno.het()] < metric[geno.homoSecond()]){
			genotypeAtMax = geno.homoSecond();
			_calledGenotype = "1/1";
			genotypeAtSecond = geno.het();
		} else {
			//het and homoSecond are equal: pick at random
			if(randomGenerator().getRand() < 0.5){
				genotypeAtMax = geno.het();
				_calledGenotype = "0/1";
				genotypeAtSecond = geno.homoSecond();
			} else {
				genotypeAtMax = geno.homoSecond();
				_calledGenotype = "1/1";
				genotypeAtSecond = geno.het();
			}
		}
	} else {
		//homoFirst and het are equal
		if(metric[geno.homoFirst()] > metric[geno.homoSecond()]){
			//pick at random between homoFirst and het
			if(randomGenerator().getRand() < 0.5){
				genotypeAtMax = geno.homoFirst();
				_calledGenotype = "0/0";
				genotypeAtSecond = geno.het();
			} else {
				genotypeAtMax = geno.het();
				_calledGenotype = "0/1";
				genotypeAtSecond = geno.homoFirst();
			}
		} else if(metric[geno.homoFirst()] < metric[geno.homoSecond()]){
			genotypeAtMax = geno.homoSecond();
			//homoFirst and het are equal: pick at random
			if(randomGenerator().getRand() < 0.5){
				genotypeAtSecond = geno.homoFirst();
			} else {
				genotypeAtSecond = geno.het();
			}
		} else {
			//all are equal: pick at random
			if(randomGenerator().getRand() < 0.333333333333333){
				genotypeAtMax = geno.homoFirst();
				_calledGenotype = "0/0";
				if(randomGenerator().getRand() < 0.5){
					genotypeAtSecond = geno.het();
				} else {
					genotypeAtSecond = geno.homoSecond();
				}
			} else if(randomGenerator().getRand() < 0.66666666666666){
				genotypeAtMax = geno.het();
				_calledGenotype = "0/1";
				if(randomGenerator().getRand() < 0.5){
					genotypeAtSecond = geno.homoFirst();
				} else {
					genotypeAtSecond = geno.homoSecond();
				}
			} else {
				genotypeAtMax = geno.homoSecond();
				_calledGenotype = "1/1";
				if(randomGenerator().getRand() < 0.5){
					genotypeAtSecond = geno.het();
				} else {
					genotypeAtSecond = geno.homoFirst();
				}
			}
		}
	}

	if(!_allowKnownAllelesCallsDifferentFromBestCall){
		//check if call matches metric of best call
		if(metric[genotypeAtMax] < *std::max_element(metric.begin(), metric.end())){
			return false;
		}
	}
	return true;
};

void TCallerDiploid::calculateImbalance(const TSite & site){
	using coretools::TBinomPValue::binomPValue;
	if(!imbalanceCalculated){
		if(!_altAlleles.empty()){
			_countAlleles(site);

			if(_altAlleles.size() == 1){
				double sum = (_alleleCounts[referenceBase] + _alleleCounts[_altAlleles[0]]);
				if(referenceBase == genometools::Base::N || sum == 0){
					AB = '.'; AI = '.';
				} else {
					AB = toString(_alleleCounts[referenceBase] / sum);
					AI = toString(binomPValue(_alleleCounts[referenceBase], _alleleCounts[_altAlleles[0]]));
				}
			} else {
				if(isHeterozygous(genotypeAtMax)){
					double sum = (double) _alleleCounts[first(genotypeAtMax)] + _alleleCounts[second(genotypeAtMax)];
					if(sum == 0){
						AB = '.'; AI = '.';
					} else {
						AB = toString(_alleleCounts[first(genotypeAtMax)] / sum);
						AI = toString(binomPValue(_alleleCounts[first(genotypeAtMax)], _alleleCounts[second(genotypeAtMax)]));
					}
				} else { // is homo -> do it against the second alternative allele
					double sum = (double) _alleleCounts[_altAlleles[0]] + _alleleCounts[_altAlleles[1]];
					if(sum == 0){
						AB = '.'; AI = '.';
					} else {
						AB = toString(_alleleCounts[_altAlleles[0]] / sum);
						AI = toString(binomPValue(_alleleCounts[_altAlleles[0]], _alleleCounts[_altAlleles[1]]));
					}
				}
			}
		}
		imbalanceCalculated = true;
	}
};

std::string TCallerDiploid::_getVCFGenotypeString_AB(const TSite & site, TGenotypeLikelihoods &){
	if(_altAlleles.empty()) return ".";

	calculateImbalance(site);
	return AB;
};

std::string TCallerDiploid::_getVCFGenotypeString_AI(const TSite & site, TGenotypeLikelihoods &){
	if(_altAlleles.empty()) return ".";

	calculateImbalance(site);
	return AI;
};

/////////////////////////////////////////////////////////
// TCallerMLE
/////////////////////////////////////////////////////////
TCallerMLE::TCallerMLE():TCallerDiploid(){
	//caller settings
	_callerName = "MLE Caller";
	_filenameExtention = "_maximumLikelihood.vcf";
	logfile().list("Will use the " + _callerName + ".");

	//parse VCF fields
	_setAcceptableFields(&_VCFGenotypeFields, "AD,GQ,GL,PL");
	printGenotypeFields("GT,DP,AD,GQ,PL"); //set default tags to print
	initializeOutput();
};

bool TCallerMLE::_callGenotype(const TSite &, TGenotypeLikelihoods & genotypeLikelihoods){
	callGenotypeFromMetric(TGenotypeProbabilities::normalize(genotypeLikelihoods));
	return true;
};

bool TCallerMLE::_callGenotypeKnownAlleles(const TSite &, TGenotypeLikelihoods & genotypeLikelihoods){
	return callGenotypeFromMetricKnownAllelesUpdateIndex(TGenotypeProbabilities::normalize(genotypeLikelihoods));
};

std::string TCallerMLE::_getVCFGenotypeString_GQ(const TSite &, TGenotypeLikelihoods & genotypeLikelihoods){
	Probability error(genotypeLikelihoods[genotypeAtSecond].get() / genotypeLikelihoods[genotypeAtMax].get());
	genometools::PhredIntProbability phred(error);
	return toString(phred);
};

std::string TCallerMLE::_getVCFGenotypeString_GL(const TSite &, TGenotypeLikelihoods & genotypeLikelihoods){
	//normalize
	TGenotypeData tmp;
	for(genometools::Genotype g = genometools::Genotype::min; g < genometools::Genotype::max; ++g){
		tmp[g] = log10(genotypeLikelihoods[g].get() / genotypeLikelihoods[genotypeAtMax].get());
	}

	//get string
	return _getPerGenotypeMetricString(tmp);
};

std::string TCallerMLE::_getVCFGenotypeString_PL(const TSite &, TGenotypeLikelihoods & genotypeLikelihoods){
	//normalize
	using genometools::Genotype;
	coretools::TStrongArray<genometools::PhredIntProbability, Genotype, 10> PL;
	genometools::PhredProbability phredMax(genotypeLikelihoods[genotypeAtMax]);
	for(genometools::Genotype g = genometools::Genotype::min; g < genometools::Genotype::max; ++g){
		PL[g] = genometools::PhredProbability(genotypeLikelihoods[g] / genotypeLikelihoods[genotypeAtMax]);
	}

	//get string
	return _getPerGenotypeMetricString(PL);
};

//------------------------------------------------------
// TCallerBayes
//------------------------------------------------------
TCallerBayes::TCallerBayes():TCallerDiploid(){
	//caller settings
	_callerName = "Bayesian Caller";
	_filenameExtention = "_maximumAPosteriori.vcf";
	_usesPrior = true;
	logfile().list("Will use the " + _callerName + ".");

	//parse VCF fields
	_setAcceptableFields(&_VCFGenotypeFields, "AD,GQ,GP");
	printGenotypeFields("GT,DP,AD,GQ,GP");  //set default tags to print
	initializeOutput();
};


bool TCallerBayes::_callGenotype(const TSite &, TGenotypeLikelihoods & genotypeLikelihoods){
	if(!_priorSet) UERROR("Can not call Bayesian genotypes: prior has not been set!");

	//calculate posterior probabilities
	_posterior = posterior(genotypeLikelihoods, *_genotypePrior);

	//call
	callGenotypeFromMetric(_posterior);
	return true;
};

bool TCallerBayes::_callGenotypeKnownAlleles(const TSite &, TGenotypeLikelihoods & genotypeLikelihoods){
	if(!_priorSet) UERROR("Can not call Bayesian genotypes: prior has not been set!");

	//calculate posterior probabilities
	_posterior = posterior(genotypeLikelihoods, *_genotypePrior);

	//call
	return callGenotypeFromMetricKnownAllelesUpdateIndex(_posterior);
};

std::string TCallerBayes::_getVCFGenotypeString_GQ(const TSite &, TGenotypeLikelihoods &){
	return (std::string) genometools::PhredIntProbability(_posterior[genotypeAtMax].complement());
};

std::string TCallerBayes::_getVCFGenotypeString_GP(const TSite &, TGenotypeLikelihoods &){
	//posterior to phred int
	coretools::TStrongArray<genometools::PhredIntProbability, genometools::Genotype, 10> tmp;
	for(genometools::Genotype g = genometools::Genotype::min; g < genometools::Genotype::max; ++g){
		tmp[g] = _posterior[g];
	}

	//get string
	return _getPerGenotypeMetricString(tmp);
};

//------------------------------------------------------
// TCall
// the class to perform calls based on windows
//------------------------------------------------------
TCall::TCall():TGenome_windows(){
	//initialize caller
	logfile().startIndent("Initializing caller:");
	std::string method = parameters().getParameterWithDefault<std::string>("method", "MLE");
	if(method == "randomBase"){
		_caller = std::make_unique<TCallerRandomBase>();
	} else if(method == "majorityBase"){
		_caller = std::make_unique<TCallerMajorityBase>();
	} else if(method == "consensify"){
		_caller = std::make_unique<TCallerConsensify>(_downsampleDepth);
	} else if(method == "allelePresence"){
		_caller = std::make_unique<TCallerAllelePresence>();
	} else if(method == "MLE"){
		_caller = std::make_unique<TCallerMLE>();
	} else if(method == "Bayesian"){
		_caller = std::make_unique<TCallerBayes>();
	} else if(method == "gVCF"){
		UERROR("GVCF NOT YET IMPLEMENTED!");
		_caller->printSitesWithNoData();
	} else UERROR("Unknown calling method '", method, "'! Use randomBase, allelePresence, MLE, Bayesian or gVCF.");
	logfile().list("Will use the " + _caller->name() + ". (parameter 'method')");

	//prior setting
	if(_caller->usesPrior()){
		_initializeGenotypePrior();
	} else {
		_prior = std::make_unique<TGenotypePriorUniform>();
	}
	_caller->setPrior(_prior->getPointerToPrior());

	//read output settings
	_caller->initializeOutput();

	//open output file
	std::string sampleName = parameters().getParameterWithDefault<std::string>("sampleName", _outputName);
	logfile().list("Will use sample name '" + sampleName + "'. (parameter 'sampleName')");
	_caller->openVCF(_outputName + "_calls", sampleName);

	//limit to sites with known alleles?
	if(parameters().parameterExists("alleles")){
		logfile().startIndent("Will limit calls to sites with known alleles (parameter 'alleles'):");
		_openSiteSubset("alleles");
		logfile().endIndent();
	} else {
		logfile().list("Will call without prior knowledge on alleles. (use 'alleles' to provide known alleles)");
		//make sure FASTA is open unless alleles are provided
		_openReference(true);
	}
	logfile().endIndent();
};

void TCall::_initializeGenotypePrior(){
	logfile().startIndent("Initializing genotype prior:");
	//read prior from parameters
	std::string priorMethod = parameters().getParameterWithDefault<std::string>("prior", "theta");
	if(priorMethod == "unif"){
		_prior = std::make_unique<TGenotypePriorUniform>();
		logfile().list("Will use a uniform prior with equal weights for all genotypes.");
	} else if(priorMethod == "theta"){
		if(parameters().parameterExists("fixedTheta")){
			double theta = parameters().getParameter<double>("fixedTheta");
			logfile().list("Will use a fixed theta = " + toString(theta));
			bool equalBaseFreq = parameters().parameterExists("equalBaseFreq");
			if(equalBaseFreq)
				logfile().list("Will use equal base frequencies.");
			else
				logfile().list("Will estimate base frequencies individually for each window.");
			_prior = std::make_unique<TGenotypePriorFixedTheta>(theta, equalBaseFreq);
		} else {
			logfile().list("Will use a prior based on theta and base frequencies estimated individually for each window.");
			std::string thetaOuputName = _outputName + "_theta_estimates.txt.gz";
			if(parameters().parameterExists("defaultTheta")){
				double defaultTheta = parameters().getParameter<double>("defaultTheta");
				logfile().list("Will use a default theta of ", defaultTheta, " for windows with limited data.");
				_prior = std::make_unique<TGenotypePriorTheta>(thetaOuputName, defaultTheta);
			} else
				_prior = std::make_unique<TGenotypePriorTheta>(thetaOuputName);
		}
	} else UERROR("Unknown prior type '", priorMethod, "'!");
	logfile().endIndent();
};

void TCall::_call(){
	uint32_t pos = 0;
	for(auto& s : _window){
		_genoLik = _genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(s);
		_caller->call(_window.chrName(), _window.positionOnChr(pos), s, _genoLik);
		++pos;
	}
};

void TCall::_callKnwonAlleles(){
	//check if we need to process this window
	if(_subset->hasPositionsInWindow(_window)){
		//add reference to sites
		_window.addReferenceBaseToSites(*_subset);

		//only run over sites listed in that window
		auto thesePositions = _subset->getPositionInWindow(_window);
		for(auto& it : thesePositions){
			//calculate genotype likelihoods
			uint32_t internalPos = it - _window.from();
			TSite& site = _window[internalPos];
			site.refBase = it.ref();
			_genoLik = _genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(site);
			_caller->call(_window.chrName(), _window.positionOnChr(internalPos), site, _genoLik, it.ref(), it.alt());
		}
	}
};


void TCall::_handleWindow(){
	if(_window.passedFilters() || _caller->printSitesWithNoData()){
		//update genotype prior
		_prior->update(_window, _genotypeLikelihoodCalculator);

		//call
		logfile().listFlushTime("Calling genotypes ...");
		if(_subset){
			_callKnwonAlleles();
		} else {
			_call();
		}
		logfile().doneTime();
	}
};

void TCall::run(){
	_traverseBAMWindows();
};

}; // end namespace
