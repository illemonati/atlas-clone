/*
 * TCaller.cpp
 *
 *  Created on: Nov 17, 2018
 *      Author: phaentu
 */

#include "TCaller.h"
#include "GenotypeFunctions.h"
#include "TSite.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Strings/fillContainer.h"
#include "coretools/Strings/toString.h"
#include "coretools/Types/probability.h"
#include "genometools/Genotypes/Base.h"
#include "genometools/Genotypes/BiallelicGenotypesWithAlleles.h"

namespace GenomeTasks{

using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using coretools::Probability;
using coretools::P;
using coretools::str::fillContainerFromStringAny;
using coretools::str::toString;
using coretools::user_assert;

using genometools::Base;
using genometools::TGenotypeLikelihoods;
using genometools::TBaseCounts;
using genometools::TGenotypeProbabilities;

using GenotypeLikelihoods::TSite;
using GenotypeLikelihoods::posterior;

namespace impl {
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
	user_assert(fields->useField(tag), "VCF ", fields->type(), " field '", tag, "' can not be printed by the ", _callerName, "!");
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
	if(parameters().exists("infoFields")){
		logfile().listFlush("Parsing VCF info fields ...");
		printInfoFields(parameters().get<std::string>("infoFields"));
		logfile().done();
	}
	if(_VCFInfoFields.numUsed() > 0){
		logfile().list("Will print these VCF info fields: " + _VCFInfoFields.getListOfUsedFields(", ") + ". (parameter 'infoFields')");
	} else {
		logfile().list("Will not print any VCF info fields. (parameter 'infoFields')");
	}

	//genotype fields
	if(parameters().exists("formatFields")){
		logfile().listFlush("Parsing VCF format fields ...");
		printGenotypeFields(parameters().get<std::string>("formatFields"));
		logfile().done();
	}
	if(_VCFGenotypeFields.numUsed() > 0){
		logfile().list("Will print these VCF format fields: " + _VCFGenotypeFields.getListOfUsedFields(", ") + ". (parameter 'formatFields')");
	} else {
		logfile().list("Will not print any VCF format fields. (parameter 'formatFields')");
	}

	//other output options
	if(parameters().exists("printAll")){
		_printSitesWithNoData = true;
		logfile().list("Will print all sites, also those without data. (parameter 'printAll')");
	} else {
		_printSitesWithNoData = false;
		logfile().list("Will print only sites with data. (use 'printAll' to print all);");
	}

	if(parameters().exists("noAltIfHomoRef")){
		_printAltIfHomoRef = false;
		logfile().list("Will not print an alternative allele if the call is homozygous reference. (parameter 'noAltIfHomoRef')");
	} else {
		_printAltIfHomoRef = true;
		logfile().list("Will print the most likely alternative allele even if the call is homozygous reference. (use 'noAltIfHomoRef' to suppress)");
	}

	if(parameters().exists("noTriallelic")){
		_allowTriallelicSites = false;
		logfile().list("Will not call genotypes resulting in two alternative alleles. (parameter 'noTriallelic')");
	} else {
		_allowTriallelicSites = true;
		logfile().list("Will allow for genotypes with two alternative alleles. (use 'noTriallelic' to suppress)");
	}

	if(parameters().exists("noCallsViolatingBest")){
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
	_vcf.open(_filename);

	//write header
	_writeVCFHeader(sampleName);
}

void TCaller::closeVCF() {
	if (_vcf.isOpen()) _vcf.close();
}

void TCaller::_writeVCFHeader(const std::string & sampleName){
	//write header
	_vcf.writeln("##fileformat=VCFv4.2");
	_vcf.writeln("##source=atlas");

	//write INFO and GENOTYPE fields
	_VCFInfoFields.writeVCFHeader(_vcf);
	_VCFGenotypeFields.writeVCFHeader(_vcf);

	//write column header
	_vcf.writeln("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",sampleName);
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
		else throw coretools::TUserError("No function defined for VCF ", _VCFInfoFields.type(), " field '", *it, "'! @Programmer: add function to TTCaller::fillInfoFieldFunctionPointers()!");
	}

};

	std::string TCaller::_getVCFInfoString_DP(const TSite & site, const genometools::TGenotypeLikelihoods &){
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
		else throw coretools::TUserError("No function defined for VCF ", _VCFGenotypeFields.type(), " field '", *it, "'! @Programmer: add function to TTCaller::fillGenotypeFieldFunctionPointers()!");

		//add to format string
		if(_genotypeFormatString.length() > 0) _genotypeFormatString += ':';
		_genotypeFormatString += *it;
	}
};

std::string TCaller::_getVCFGenotypeString_GT(const TSite &, const TGenotypeLikelihoods &){
	return _calledGenotype;
};

std::string TCaller::_getVCFGenotypeString_DP(const TSite & site, const TGenotypeLikelihoods &){
	return toString(site.depth());
};

std::string TCaller::_getVCFGenotypeString_AD(const TSite & site, const TGenotypeLikelihoods &){
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
std::string TCaller::_composeVCFString(std::vector<std::string (TCaller::*)(const TSite & site, const TGenotypeLikelihoods & genotypeLikelihoods)> & vec, const TSite & site, const TGenotypeLikelihoods & genotypeLikelihoods){
	//no info fields?
	if(vec.empty()) return ".";

	//add first info
	auto it = vec.begin();
	std::string info = (this->*(*it))(site, genotypeLikelihoods);
	++it;

	//loop over rest
	for(; it != vec.end(); ++it) {
		info += ':' + (this->*(*it))(site, genotypeLikelihoods);
	}

	return info;
};

void TCaller::_writeAlternativeAllelesToVCF(){
	if(_altAlleles.size() == 0){
		_vcf.write(".");
	} else {
		_vcf.write(_altAlleles[0]);
		for(size_t i=1; i<_altAlleles.size(); ++i)
			_vcf.write(',',_altAlleles[i]);
	}
};

void TCaller::_writeCallToVCF(const std::string & chr, const long pos, const TSite & site, const TGenotypeLikelihoods & genotypeLikelihoods){
	//apply filter on alternative alleles
	if(!_printAltIfHomoRef && (_calledGenotype == "0/0" || _calledGenotype == "0"))
		_altAlleles.clear();

	//write chr, position and (no) variant ID
	_vcf.write(chr,'\t',pos + 1,"\t.\t"); //all internal positions are zero-based!

	//write reference and alternative alleles
	_vcf.write(site.refBase,"\t");
	_writeAlternativeAllelesToVCF();

	//write (no) variant quality and (no) filter
	_vcf.write("\t.\t.");

	//write info fields
	_vcf.write('\t',_composeVCFString(_VCFInfoFunctionsVec, site, genotypeLikelihoods));

	//write genotype fields
	_vcf.write('\t',_genotypeFormatString,'\t',_composeVCFString(_VCFGenotypeFunctionsVec, site, genotypeLikelihoods));

	//end with new line
	_vcf.endln();

	//clean up storage
	_clearAfterCall();
};

void TCaller::_writeMissingDataToVCF(const TSite & site){
	if(_printSitesWithNoData)
		_vcf.write("\t.\t",site.refBase,"\t.\t.\t.\t.\tGT:DP\t",_missingGenotype,":0");
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
		_alleleCounts   = site.countAlleles();
		_allelesCounted = true;
	}
};

bool TCaller::_callGenotype(const TSite &, const TGenotypeLikelihoods &){
	_calledGenotype = "./.";
	return true;
};

bool TCaller::_callGenotypeKnownAlleles(const TSite &, const TGenotypeLikelihoods &){
	_calledGenotype = "./.";
	return true;
};

void TCaller::call(const std::string & chr, long pos, const TSite & site, const TGenotypeLikelihoods & genotypeLikelihoods){
	//set reference base from site
	referenceBase = site.refBase;

	//check if there is data
	if(site.empty() || (referenceBase == genometools::Base::N && !_allowTriallelicSites) || !_callGenotype(site, genotypeLikelihoods))  {
		_writeMissingDataToVCF(site);
	} else {
		_writeCallToVCF(chr, pos, site, genotypeLikelihoods);
	}
};

void TCaller::call(const std::string & chr, long pos, const TSite & site, const TGenotypeLikelihoods & genotypeLikelihoods, const genometools::Base & firstAllele, const genometools::Base & secondAllele){
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

bool TCallerRandomBase::_callGenotype(const TSite & site, const TGenotypeLikelihoods &){
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

bool TCallerRandomBase::_callGenotypeKnownAlleles(const TSite & site, const TGenotypeLikelihoods &){
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

bool TCallerMajorityBase::_callGenotype(const TSite & site, const TGenotypeLikelihoods &){
	//get per allele counts
	_countAlleles(site);

	//call majority
	const auto [majorityBase, second] = impl::sampleFirstSecond(_alleleCounts);

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

bool TCallerMajorityBase::_callGenotypeKnownAlleles(const TSite & site, const TGenotypeLikelihoods &){
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
		const auto majorityBase = randomGenerator().sampleIndexOfMaxima<genometools::TBaseCounts, Base, 4>(_alleleCounts);

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
	user_assert(DownsampleDepth > 0, "Consensify caller requires downsampling of reads! Use 'downsample' to specify.");

	_downsampleDepth = DownsampleDepth;
	_minMajorityDepth = ceil(_downsampleDepth / 2.0);

	logfile().list("Will call based on a majority out of " + toString(_downsampleDepth) + " bases. (parameter 'downsample')");
};

bool TCallerConsensify::_callGenotype(const TSite & site, const TGenotypeLikelihoods &){
	DEV_ASSERT(site.depth() <= _downsampleDepth);

	//get per allele counts
	_countAlleles(site);

	//call majority
	const auto [majorityBase, second] = impl::sampleFirstSecond(_alleleCounts);

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

bool TCallerConsensify::_callGenotypeKnownAlleles(const TSite & site, const TGenotypeLikelihoods &){
	DEV_ASSERT(site.depth() <= _downsampleDepth);

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

void TCallerAllelePresence::_fillPosteriors(const TGenotypeLikelihoods & genotypeLikelihoods){
	using genometools::Base;
	using genometools::Genotype;
	//calculate posterior probabilities
	_posterior = posterior(genotypeLikelihoods, *_genotypePrior);

	//sum for each base
	_allelePostProb[Base::A] = P(_posterior[Genotype::AA] + _posterior[Genotype::AC] + _posterior[Genotype::AG] + _posterior[Genotype::AT]);
	_allelePostProb[Base::C] = P(_posterior[Genotype::AC] + _posterior[Genotype::CC] + _posterior[Genotype::CG] + _posterior[Genotype::CT]);
	_allelePostProb[Base::G] = P(_posterior[Genotype::AG] + _posterior[Genotype::CG] + _posterior[Genotype::GG] + _posterior[Genotype::GT]);
	_allelePostProb[Base::T] = P(_posterior[Genotype::AT] + _posterior[Genotype::CT] + _posterior[Genotype::GT] + _posterior[Genotype::TT]);
};

bool TCallerAllelePresence::_callGenotype(const TSite &, const TGenotypeLikelihoods & genotypeLikelihoods){
	user_assert(_priorSet, "Can not call AllelePresence genotypes: prior has not been set!");

	//fill posteriors for each allele
	_fillPosteriors(genotypeLikelihoods);

	//find MAP
	const auto [majorityBase, second] = impl::sampleFirstSecond(_alleleCounts);
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

bool TCallerAllelePresence::_callGenotypeKnownAlleles(const TSite &, const TGenotypeLikelihoods & genotypeLikelihoods){
	user_assert(_priorSet, "Can not call AllelePresence genotypes: prior has not been set!");

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

std::string TCallerAllelePresence::_getVCFGenotypeString_GQ(const TSite &, const TGenotypeLikelihoods &){
	return (std::string) coretools::PhredInt(_allelePostProb[_MAP].complement());
};

std::string TCallerAllelePresence::_getVCFGenotypeString_AP(const TSite &, const TGenotypeLikelihoods &){
	using genometools::Base;
	std::string ret = (std::string) coretools::PhredInt(_allelePostProb[Base::A]);
	ret += ',' + (std::string) coretools::PhredInt(_allelePostProb[Base::C]);
	ret += ',' + (std::string) coretools::PhredInt(_allelePostProb[Base::G]);
	ret += ',' + (std::string) coretools::PhredInt(_allelePostProb[Base::T]);
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
	std::tie(_genotypeAtMax, _genotypeAtSecond) = impl::sampleFirstSecond(metric);

	//decide on alternative alleles
	if(first(_genotypeAtMax) == referenceBase){
		if(second(_genotypeAtMax) == referenceBase){
			_calledGenotype = "0/0";
			//MLE is homozygous reference -> find second best allele
			if(first(_genotypeAtSecond) == referenceBase)
				_altAlleles.push_back(second(_genotypeAtSecond));
			else if(second(_genotypeAtSecond) == referenceBase)
				_altAlleles.push_back(first(_genotypeAtSecond));
			else {
				//none is ref, pick at random
				const auto b = randomGenerator().getRand() < 0.5 ? first(_genotypeAtSecond) : second(_genotypeAtSecond);
				_altAlleles.push_back( b );
			}
		} else {
			_altAlleles.push_back(second(_genotypeAtMax));
			_calledGenotype = "0/1";
		}
	} else {
		if(second(_genotypeAtMax) == referenceBase){
			_altAlleles.push_back(first(_genotypeAtMax));
			_calledGenotype = "0/1";
		} else {
			if(isHomozygous(_genotypeAtMax)){
				_altAlleles.push_back(first(_genotypeAtMax));
				_calledGenotype = "1/1";

				//find second best allele, but give preference to reference in case likelihoods are equal
				if(_allowTriallelicSites){
					//int hetRef = genoMap.getGenotype(referenceBase, genoMap.genotypeToBase[indexOfMax][0]);
					//int homRef = genoMap.getGenotype(referenceBase, referenceBase);

					//only use second alternative allele in case het genotype with reference is less likely
					if(referenceBase == Base::N ||
					   (metric[genotype(referenceBase, first(_genotypeAtMax))] < metric[ _genotypeAtSecond]
						&& metric[genotype(referenceBase, referenceBase) ] < metric[ _genotypeAtSecond ])){
						if(first(_genotypeAtSecond) == referenceBase || first(_genotypeAtSecond) == _altAlleles[0])
							_altAlleles.push_back(second(_genotypeAtSecond));
						else if(second(_genotypeAtSecond) == referenceBase || second(_genotypeAtSecond) == _altAlleles[0])
							_altAlleles.push_back(first(_genotypeAtSecond));
						else {
							//none is ref, pick at random
							const auto b = randomGenerator().getRand() < 0.5 ? first(_genotypeAtSecond) : second(_genotypeAtSecond);
							_altAlleles.push_back(b);
						}
					}
				}
			} else {
				if(_allowTriallelicSites){
					//allow triallelic sites
					_altAlleles.push_back(first(_genotypeAtMax));
					_altAlleles.push_back(second(_genotypeAtMax));
					_calledGenotype = "1/2";
				} else {
					//decide on which of the two alternative alleles to pick -> check second highest
					if(isHomozygous(_genotypeAtSecond)){
						if(first(_genotypeAtSecond) == first(_genotypeAtMax))
							_altAlleles.push_back(first(_genotypeAtMax));
						else if(first(_genotypeAtSecond) == second(_genotypeAtMax))
							_altAlleles.push_back(second(_genotypeAtMax));
						else {
							//neither alt 0 nor 1, pick at random
							const auto b = randomGenerator().getRand() < 0.5 ? first(_genotypeAtMax) : second(_genotypeAtMax);
							_altAlleles.push_back(b);
						}
					} else {
						if(first(_genotypeAtSecond) == referenceBase && (second(_genotypeAtSecond) == first(_genotypeAtMax) || second(_genotypeAtSecond) == first(_genotypeAtMax))){
							_altAlleles.push_back(second(_genotypeAtSecond));
						} else if(second(_genotypeAtSecond) == referenceBase && (first(_genotypeAtSecond) == first(_genotypeAtMax) || first(_genotypeAtSecond) == second(_genotypeAtMax))){
							_altAlleles.push_back(first(_genotypeAtSecond));
						} else {
							//neither alt 0 nor 1, pick at random
							const auto b = randomGenerator().getRand() < 0.5 ? first(_genotypeAtMax) : second(_genotypeAtMax);
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
			_genotypeAtMax = geno.homoFirst();
			_calledGenotype = "0/0";
			if(metric[geno.het()] > metric[geno.homoSecond()]){
				_genotypeAtSecond = geno.het();
			} else if(metric[geno.het()] < metric[geno.homoSecond()]){
				_genotypeAtSecond = geno.homoSecond();
			} else {
				//het and homoSecond are equal: pick at random
				if(randomGenerator().getRand() < 0.5){
					_genotypeAtSecond = geno.het();
				} else {
					_genotypeAtSecond = geno.homoSecond();
				}
			}
		} else if(metric[geno.homoFirst()] < metric[geno.homoSecond()]){
			_genotypeAtMax = geno.homoSecond();
			_calledGenotype = "1/1";
			_genotypeAtSecond = geno.homoFirst();
		} else {
			//homoFirst and homoSecond are equal: pick at random
			if(randomGenerator().getRand() < 0.5){
				_genotypeAtMax = geno.homoFirst();
				_calledGenotype = "0/0";
				_genotypeAtSecond = geno.homoSecond();
			} else {
				_genotypeAtMax = geno.homoSecond();
				_calledGenotype = "1/1";
				_genotypeAtSecond = geno.homoFirst();
			}
		}
	} else if(metric[geno.homoFirst()] < metric[geno.het()]){
		if(metric[geno.het()] > metric[geno.homoSecond()]){
			_genotypeAtMax = geno.het();
			_calledGenotype = "0/1";
			if(metric[geno.homoFirst()] > metric[geno.homoSecond()]){
				_genotypeAtSecond = geno.homoFirst();
			} else if(metric[geno.homoFirst()] < metric[geno.homoSecond()]){
				_genotypeAtSecond = geno.homoSecond();
			} else {
				//homoFirst and homoSecond are equal: pick at random
				if(randomGenerator().getRand() < 0.5){
					_genotypeAtSecond = geno.homoFirst();
				} else {
					_genotypeAtSecond = geno.homoSecond();
				}
			}
		} else if(metric[geno.het()] < metric[geno.homoSecond()]){
			_genotypeAtMax = geno.homoSecond();
			_calledGenotype = "1/1";
			_genotypeAtSecond = geno.het();
		} else {
			//het and homoSecond are equal: pick at random
			if(randomGenerator().getRand() < 0.5){
				_genotypeAtMax = geno.het();
				_calledGenotype = "0/1";
				_genotypeAtSecond = geno.homoSecond();
			} else {
				_genotypeAtMax = geno.homoSecond();
				_calledGenotype = "1/1";
				_genotypeAtSecond = geno.het();
			}
		}
	} else {
		//homoFirst and het are equal
		if(metric[geno.homoFirst()] > metric[geno.homoSecond()]){
			//pick at random between homoFirst and het
			if(randomGenerator().getRand() < 0.5){
				_genotypeAtMax = geno.homoFirst();
				_calledGenotype = "0/0";
				_genotypeAtSecond = geno.het();
			} else {
				_genotypeAtMax = geno.het();
				_calledGenotype = "0/1";
				_genotypeAtSecond = geno.homoFirst();
			}
		} else if(metric[geno.homoFirst()] < metric[geno.homoSecond()]){
			_genotypeAtMax = geno.homoSecond();
			//homoFirst and het are equal: pick at random
			if(randomGenerator().getRand() < 0.5){
				_genotypeAtSecond = geno.homoFirst();
			} else {
				_genotypeAtSecond = geno.het();
			}
		} else {
			//all are equal: pick at random
			if(randomGenerator().getRand() < 0.333333333333333){
				_genotypeAtMax = geno.homoFirst();
				_calledGenotype = "0/0";
				if(randomGenerator().getRand() < 0.5){
					_genotypeAtSecond = geno.het();
				} else {
					_genotypeAtSecond = geno.homoSecond();
				}
			} else if(randomGenerator().getRand() < 0.66666666666666){
				_genotypeAtMax = geno.het();
				_calledGenotype = "0/1";
				if(randomGenerator().getRand() < 0.5){
					_genotypeAtSecond = geno.homoFirst();
				} else {
					_genotypeAtSecond = geno.homoSecond();
				}
			} else {
				_genotypeAtMax = geno.homoSecond();
				_calledGenotype = "1/1";
				if(randomGenerator().getRand() < 0.5){
					_genotypeAtSecond = geno.het();
				} else {
					_genotypeAtSecond = geno.homoFirst();
				}
			}
		}
	}

	if(!_allowKnownAllelesCallsDifferentFromBestCall){
		//check if call matches metric of best call
		if(metric[_genotypeAtMax] < *std::max_element(metric.begin(), metric.end())){
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
				if(isHeterozygous(_genotypeAtMax)){
					double sum = (double) _alleleCounts[first(_genotypeAtMax)] + _alleleCounts[second(_genotypeAtMax)];
					if(sum == 0){
						AB = '.'; AI = '.';
					} else {
						AB = toString(_alleleCounts[first(_genotypeAtMax)] / sum);
						AI = toString(binomPValue(_alleleCounts[first(_genotypeAtMax)], _alleleCounts[second(_genotypeAtMax)]));
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

std::string TCallerDiploid::_getVCFGenotypeString_AB(const TSite & site, const TGenotypeLikelihoods &){
	if(_altAlleles.empty()) return ".";

	calculateImbalance(site);
	return AB;
};

std::string TCallerDiploid::_getVCFGenotypeString_AI(const TSite & site, const TGenotypeLikelihoods &){
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

bool TCallerMLE::_callGenotype(const TSite &, const TGenotypeLikelihoods & genotypeLikelihoods){
	callGenotypeFromMetric(TGenotypeProbabilities::normalize(genotypeLikelihoods));
	return true;
};

bool TCallerMLE::_callGenotypeKnownAlleles(const TSite &, const TGenotypeLikelihoods & genotypeLikelihoods){
	return callGenotypeFromMetricKnownAllelesUpdateIndex(TGenotypeProbabilities::normalize(genotypeLikelihoods));
};

std::string TCallerMLE::_getVCFGenotypeString_GQ(const TSite &, const TGenotypeLikelihoods & genotypeLikelihoods){
	Probability error(genotypeLikelihoods[_genotypeAtSecond].get() / genotypeLikelihoods[_genotypeAtMax].get());
	coretools::PhredInt phred(error);
	return toString(phred);
};

std::string TCallerMLE::_getVCFGenotypeString_GL(const TSite &, const TGenotypeLikelihoods & genotypeLikelihoods){
	//normalize
	genometools::TGenotypeData tmp;
	for(genometools::Genotype g = genometools::Genotype::min; g < genometools::Genotype::max; ++g){
		tmp[g] = log10(genotypeLikelihoods[g].get() / genotypeLikelihoods[_genotypeAtMax].get());
	}

	//get string
	return _getPerGenotypeMetricString(tmp);
};

std::string TCallerMLE::_getVCFGenotypeString_PL(const TSite &, const TGenotypeLikelihoods & genotypeLikelihoods){
	//normalize
	using genometools::Genotype;
	coretools::TStrongArray<coretools::PhredInt, Genotype, 10> PL;
	for(genometools::Genotype g = genometools::Genotype::min; g < genometools::Genotype::max; ++g){
		auto p = genotypeLikelihoods[g];
		p.scale(genotypeLikelihoods[_genotypeAtMax]);
		PL[g] = coretools::PhredInt(p);
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


bool TCallerBayes::_callGenotype(const TSite &, const TGenotypeLikelihoods & genotypeLikelihoods){
	user_assert(_priorSet, "Can not call Bayesian genotypes: prior has not been set!");

	//calculate posterior probabilities
	_posterior = posterior(genotypeLikelihoods, *_genotypePrior);

	//call
	callGenotypeFromMetric(_posterior);
	return true;
};

bool TCallerBayes::_callGenotypeKnownAlleles(const TSite &, const TGenotypeLikelihoods & genotypeLikelihoods){
	user_assert(_priorSet, "Can not call Bayesian genotypes: prior has not been set!");

	//calculate posterior probabilities
	_posterior = posterior(genotypeLikelihoods, *_genotypePrior);

	//call
	return callGenotypeFromMetricKnownAllelesUpdateIndex(_posterior);
};

std::string TCallerBayes::_getVCFGenotypeString_GQ(const TSite &, const TGenotypeLikelihoods &){
	return (std::string) coretools::PhredInt(_posterior[_genotypeAtMax].complement());
};

std::string TCallerBayes::_getVCFGenotypeString_GP(const TSite &, const TGenotypeLikelihoods &){
	//posterior to phred int
	coretools::TStrongArray<coretools::PhredInt, genometools::Genotype, 10> tmp;
	for(genometools::Genotype g = genometools::Genotype::min; g < genometools::Genotype::max; ++g){
		tmp[g] = coretools::PhredInt(_posterior[g]);
	}

	//get string
	return _getPerGenotypeMetricString(tmp);
};

//------------------------------------------------------
// TCall
// the class to perform calls based on windows
//------------------------------------------------------
TCall::TCall() {
	//initialize caller
	logfile().startIndent("Initializing caller:");
	std::string method = parameters().get<std::string>("method", "MLE");
	if(method == "randomBase"){
		_caller = std::make_unique<TCallerRandomBase>();
	} else if(method == "majorityBase"){
		_caller = std::make_unique<TCallerMajorityBase>();
	} else if(method == "consensify"){
		const auto dDepth = parameters().get<size_t>("downsample");
		_caller = std::make_unique<TCallerConsensify>(dDepth);
	} else if(method == "allelePresence"){
		_caller = std::make_unique<TCallerAllelePresence>();
	} else if(method == "MLE"){
		_caller = std::make_unique<TCallerMLE>();
	} else if(method == "Bayesian"){
		_caller = std::make_unique<TCallerBayes>();
	} else if(method == "gVCF"){
		throw coretools::TUserError("GVCF NOT YET IMPLEMENTED!");
		_caller->printSitesWithNoData();
	} else throw coretools::TUserError("Unknown calling method '", method, "'! Use randomBase, allelePresence, MLE, Bayesian or gVCF.");
	logfile().list("Will use the " + _caller->name() + ". (parameter 'method')");

	//prior setting
	if(_caller->usesPrior()){
		_initializeGenotypePrior();
	} else {
		_prior = std::make_unique<GenotypeLikelihoods::TGenotypePriorUniform>();
	}
	_caller->setPrior(_prior->getPointerToPrior());

	//read output settings
	_caller->initializeOutput();

	//open output file
	std::string sampleName = parameters().get<std::string>("sampleName", _genome.outputName());
	logfile().list("Will use sample name '" + sampleName + "'. (parameter 'sampleName')");
	_caller->openVCF(_genome.outputName() + "_calls", sampleName);

	//limit to sites with known alleles?
	if(parameters().exists("alleles")){
		_windows.openSiteSubset("alleles", _genome.bamFile().chromosomes(), genometools::Morphic::Both);
		_RefN = parameters().exists("allowRefN") ? 0 : -1;
		if (_RefN >= 0) {
			logfile().list("Allowing cases where Ref = N and Alt != N, treated as N/N. (parameter 'allowRefN')");
		} else {
			logfile().list("Not allowing cases where Ref = N and Alt != N. (use 'allowRefN' to allow)");
		}
	} else {
		logfile().list("Will call without prior knowledge on alleles. (use 'alleles' to provide known alleles)");
		//make sure FASTA is open unless alleles are provided
		_windows.requireReference();
	}
	logfile().endIndent();
};

void TCall::_initializeGenotypePrior(){
	logfile().startIndent("Initializing genotype prior:");
	//read prior from parameters
	std::string priorMethod = parameters().get<std::string>("prior", "theta");
	if(priorMethod == "unif"){
		_prior = std::make_unique<GenotypeLikelihoods::TGenotypePriorUniform>();
		logfile().list("Will use a uniform prior with equal weights for all genotypes.");
	} else if(priorMethod == "theta"){
		if(parameters().exists("fixedTheta")){
			double theta = parameters().get<double>("fixedTheta");
			logfile().list("Will use a fixed theta = " + toString(theta));
			bool equalBaseFreq = parameters().exists("equalBaseFreq");
			if(equalBaseFreq)
				logfile().list("Will use equal base frequencies.");
			else
				logfile().list("Will estimate base frequencies individually for each window.");
			_prior = std::make_unique<GenotypeLikelihoods::TGenotypePriorFixedTheta>(theta, equalBaseFreq);
		} else {
			logfile().list("Will use a prior based on theta and base frequencies estimated individually for each window.");
			std::string thetaOuputName = _genome.outputName() + "_theta_estimates.txt.gz";
			if(parameters().exists("defaultTheta")){
				double defaultTheta = parameters().get<double>("defaultTheta");
				logfile().list("Will use a default theta of ", defaultTheta, " for windows with limited data.");
				_prior = std::make_unique<GenotypeLikelihoods::TGenotypePriorTheta>(thetaOuputName, defaultTheta);
			} else
				_prior = std::make_unique<GenotypeLikelihoods::TGenotypePriorTheta>(thetaOuputName);
		}
	} else throw coretools::TUserError("Unknown prior type '", priorMethod, "'!");
	logfile().endIndent();
};

void TCall::_call(GenotypeLikelihoods::TWindow& window){
	uint32_t pos = 0;
	for(auto& s : window){
		const auto genoLik = _genome.errorModels().calculateGenotypeLikelihoods(s);
		_caller->call(window.chrName(), window.positionOnChr(pos), s, genoLik);
		++pos;
	}
};

void TCall::_callKnwonAlleles(GenotypeLikelihoods::TWindow& window){
	//check if we need to process this window
	const auto& alleles = _windows.alleles();
	if(alleles.overlaps(window)){
		//add reference to sites
		window.addReferenceBaseToSites(alleles);

		//only run over sites listed in that window
		for (auto it = alleles.begin(window); it != alleles.end() && window.within(it->position); ++it) {
			//calculate genotype likelihoods
			uint32_t internalPos = it->position - window.from();
			TSite& site = window[internalPos];
			site.refBase = it->ref;
			const auto genoLik = _genome.errorModels().calculateGenotypeLikelihoods(site);

			if (it->ref == Base::N && it->alt != Base::N) {
				if (_RefN >= 0) {
					++_RefN;
				} else {
					throw coretools::TUserError("In position ", it->position.asFormattedString(_genome.bamFile().chromosomes()),
						   ", Ref = N but Alt = ", it->alt,
						   ", this case is not supported! Use parameter 'allowRefN' to treat as N/N");
				}
			}

			if (it->alt == Base::N || it->ref == Base::N) _caller->call(window.chrName(), window.positionOnChr(internalPos), site, genoLik);
			else _caller->call(window.chrName(), window.positionOnChr(internalPos), site, genoLik, it->ref, it->alt);
		}
	}
};

void TCall::_handleWindow(GenotypeLikelihoods::TWindow& window){
	if(window.passedFilters() || _caller->printSitesWithNoData()){
		//update genotype prior
		_prior->update(window, _genome.errorModels());

		//call
		logfile().listFlushTime("Calling genotypes ...");
		if(_windows.alleles()){
			_callKnwonAlleles(window);
		} else {
			_call(window);
		}
		logfile().doneTime();
	}
};

void TCall::run(){
	_traverseBAMWindows();
	if (_RefN >= 0) {
		logfile().list("Number of sites where Ref = N and Alt != N: ", _RefN, ".");
	}
};

}; // end namespace
