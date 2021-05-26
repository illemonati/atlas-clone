/*
 * TCaller.cpp
 *
 *  Created on: Nov 17, 2018
 *      Author: phaentu
 */


#include "TCaller.h"

namespace GenomeTasks{

using namespace GenotypeLikelihoods;

/////////////////////////////////////////////////////////
// TCaller
/////////////////////////////////////////////////////////
TCaller::TCaller(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator){
	_callerName = "default caller";
	_filenameExtention = ".vcf.gz";
	_randomGenerator = RandomGenerator;

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
void TCaller::_setAcceptableFields(VCF::TVCFFieldVector* fields, std::string tags){
	std::vector<std::string> vec;
	fillContainerFromStringAny(tags, vec, ",", true);
	for(std::vector<std::string>::iterator it = vec.begin(); it != vec.end(); ++it)
		fields->acceptField(*it);
};

void TCaller::_printField(VCF::TVCFFieldVector* fields, std::string tag){
	if(!fields->useField(tag))
		throw "VCF " +  fields->type() + " field '" + tag + "' can not be printed by the " + _callerName +"!";
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

void TCaller::initializeOutput(TParameters & Parameters, TLog* Logfile){
	//info fields
	if(Parameters.parameterExists("infoFields")){
		Logfile->listFlush("Parsing VCF info fields ...");
		printInfoFields(Parameters.getParameter<std::string>("infoFields"));
		Logfile->done();
	}
	if(_VCFInfoFields.numUsed() > 0){
		Logfile->list("Will print these VCF info fields: " + _VCFInfoFields.getListOfUsedFields(", ") + ". (parameter 'infoFields')");
	} else {
		Logfile->list("Will not print any VCF info fields. (parameter 'infoFields')");
	}

	//genotype fields
	if(Parameters.parameterExists("formatFields")){
		Logfile->listFlush("Parsing VCF format fields ...");
		printGenotypeFields(Parameters.getParameter<std::string>("formatFields"));
		Logfile->done();
	}
	if(_VCFGenotypeFields.numUsed() > 0){
		Logfile->list("Will print these VCF format fields: " + _VCFGenotypeFields.getListOfUsedFields(", ") + ". (parameter 'formatFields')");
	} else {
		Logfile->list("Will not print any VCF format fields. (parameter 'formatFields')");
	}

	//other output options
	if(Parameters.parameterExists("printAll")){
		_printSitesWithNoData = true;
		Logfile->list("Will print all sites, also those without data. (parameter 'printAll')");
	} else {
		_printSitesWithNoData = false;
		Logfile->list("Will print only sites with data. (use 'printAll' to print all);");
	}

	if(Parameters.parameterExists("noAltIfHomoRef")){
		_printAltIfHomoRef = false;
		Logfile->list("Will not print an alternative allele if the call is homozygous reference. (parameter 'noAltIfHomoRef')");
	} else {
		_printAltIfHomoRef = true;
		Logfile->list("Will print the most likely alternative allele even if the call is homozygous reference. (use 'noAltIfHomoRef' to suppress)");
	}

	if(Parameters.parameterExists("noTriallelic")){
		_allowTriallelicSites = false;
		Logfile->list("Will not call genotypes resulting in two alternative alleles. (parameter 'noTriallelic')");
	} else {
		_allowTriallelicSites = true;
		Logfile->list("Will allow for genotypes with two alternative alleles. (use 'noTriallelic' to suppress)");
	}

	if(Parameters.parameterExists("noCallsViolatingBest")){
		_allowKnownAllelesCallsDifferentFromBestCall = false;
		Logfile->list("Will not call genotypes from known alleles that conflict with best call across all genotypes. (parameter 'noCallsViolatingBest')");
	} else {
		_allowKnownAllelesCallsDifferentFromBestCall = true;
		Logfile->list("Will call genotypes from known alleles even if they differ from best call across all genotypes. (use 'noCallsViolatingBest' to suppress)");
	}
}

//-------------------------------------------------------------------------------------------
// open / close VCF file, print header
//-------------------------------------------------------------------------------------------
void TCaller::openVCF(const std::string FilenameTag, const std::string sampleName, TLog* logfile){
	_filename = FilenameTag  + _filenameExtention + ".gz";
	logfile->list("Writing calls to VCF file '" + _filename + "'.");
	_vcf.open(_filename.c_str());
	if(!_vcf) throw "Failed to open VCF file '" + _filename + "' for writing!";
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
		else throw "No function defined for VCF " + _VCFInfoFields.type() + " field '" + *it + "'! @Programmer: add function to TTCaller::fillInfoFieldFunctionPointers()!";
	}

};

std::string TCaller::_getVCFInfoString_DP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
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
		else throw "No function defined for VCF " + _VCFGenotypeFields.type() + " field '" + *it + "'! @Programmer: add function to TTCaller::fillGenotypeFieldFunctionPointers()!";

		//add to format string
		if(_genotypeFormatString.length() > 0) _genotypeFormatString += ':';
		_genotypeFormatString += *it;
	}
};

std::string TCaller::_getVCFGenotypeString_GT(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	return _calledGenotype;
};

std::string TCaller::_getVCFGenotypeString_DP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	return toString(site.depth());
};

std::string TCaller::_getVCFGenotypeString_AD(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	_countAlleles(site);
	std::string ret;
	if(referenceBase == BAM::N) ret = "0";
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
	_vcf << site.refBase() << "\t";
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
		_vcf << "\t.\t" << site.refBase() << "\t.\t.\t.\t.\tGT:DP\t" << _missingGenotype << ":0";
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

bool TCaller::_callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	_calledGenotype = "./.";
	return true;
};

bool TCaller::_callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	_calledGenotype = "./.";
	return true;
};

void TCaller::call(const std::string & chr, const long & pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//set reference base from site
	referenceBase = site.refBase();

	//check if there is data
	if(site.empty() || (referenceBase == BAM::N && !_allowTriallelicSites) || !_callGenotype(site, genotypeLikelihoods))
		_writeMissingDataToVCF(site);
	else {
		_writeCallToVCF(chr, pos, site, genotypeLikelihoods);
	}
};

void TCaller::call(const std::string & chr, const long & pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods, const BAM::Base & firstAllele, const BAM::Base & secondAllele){
	//check if there is data
	if(!site.empty()){
		//set reference base from site
		referenceBase = site.refBase();

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
TCallerRandomBase::TCallerRandomBase(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TCaller(Parameters, Logfile, RandomGenerator){
	//caller settings
	_callerName = "Random Base Caller";
	_filenameExtention = "_randomBase.vcf";
	Logfile->list("Will use the " + _callerName + ".");

	//parse VCF fields
	_setAcceptableFields(&_VCFInfoFields, "DP");
	_setAcceptableFields(&_VCFGenotypeFields, "GT,DP,AD");
	initializeOutput(Parameters, Logfile);

};

bool TCallerRandomBase::_callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//randomly pick a base
	BAM::Base allele = site[_randomGenerator->sample(site.depth())].base;

	//decide on alt
	if(allele == referenceBase){
		_calledGenotype = "0";
	} else {
		_altAlleles.push_back(allele);
		_calledGenotype = "1";
	}

	return true;
};

bool TCallerRandomBase::_callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(_allowKnownAllelesCallsDifferentFromBestCall){
		//pick among known alleles only
		_countAlleles(site);
		double probRef = (double) _alleleCounts[referenceBase] / (double) (_alleleCounts[referenceBase] + _alleleCounts[_altAlleles[0]]);

		//pick among known alleles
		if(_randomGenerator->getRand() < probRef){
			_calledGenotype = "0";
		} else {
			_calledGenotype = "1";
		}
		return true;
	} else {
		//pick among all alleles and check
		BAM::Base allele = site[_randomGenerator->sample(site.depth())].base;
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
TCallerMajorityBase::TCallerMajorityBase(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TCaller(Parameters, Logfile, RandomGenerator){
	//general caller settings
	_callerName = "Majority Base Caller";
	_filenameExtention = "_majorityBase.vcf";
	Logfile->list("Will use the " + _callerName + ".");

	//parse VCF fields
	initializeOutput(Parameters, Logfile);
};

bool TCallerMajorityBase::_callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//get per allele counts
	_countAlleles(site);

	//call majority
	BAM::Base majorityBase = static_cast<BAM::BaseEnum>( _pickIndexWithHighestMetric(_alleleCounts) );

	//decide on alt
	if(majorityBase == referenceBase){
		_calledGenotype = "0";

		//find second most common as alternative allele
		BAM::Base second = static_cast<BAM::BaseEnum>( _pickIndexWithSecondHighestMetric(_alleleCounts, majorityBase.get()) );
		_altAlleles.push_back(second);
	} else {
		_altAlleles.push_back(majorityBase);
		_calledGenotype = "1";
	}

	return true;
};

bool TCallerMajorityBase::_callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
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
			if(_randomGenerator->getRand() < 0.5){
				_calledGenotype = "0";
			} else {
				_calledGenotype = "1";
			}
		}
		return true;
	} else {
		//pick among all alleles and check
		//call majority
		BAM::Base majorityBase = static_cast<BAM::BaseEnum>( _pickIndexWithHighestMetric(_alleleCounts) );

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
TCallerConsensify::TCallerConsensify(const uint32_t & DownsampleDepth, TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TCaller(Parameters, Logfile, RandomGenerator){
	//general caller settings
	_callerName = "Consensify Base Caller";
	_filenameExtention = "_consensifyBase.vcf";
	Logfile->list("Will use the " + _callerName + ".");

	//parse VCF fields
	initializeOutput(Parameters, Logfile);

	//check downsample depth
	if(DownsampleDepth < 1){
		throw "Consensify caller requires downsampling of reads! Use 'downsample' to specify.";
	}
	_downsampleDepth = DownsampleDepth;
	_minMajorityDepth = ceil(_downsampleDepth / 2.0);

	Logfile->list("Will call based on a majority out of " + toString(_downsampleDepth) + " bases. (parameter 'downsample')");
};

bool TCallerConsensify::_callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(site.depth() > _downsampleDepth){
		throw std::runtime_error("bool TCallerConsensify::_callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods): depth > _downsampleDepth!");
	}

	//get per allele counts
	_countAlleles(site);

	//call majority
	BAM::Base majorityBase = static_cast<BAM::BaseEnum>( _pickIndexWithHighestMetric(_alleleCounts) );

	//check if we have sufficient depth to call
	if(_alleleCounts[majorityBase] < _minMajorityDepth){
		return false;
	}

	//decide on alt
	if(majorityBase == referenceBase){
		_calledGenotype = "0";

		//find second most common as alternative allele
		BAM::Base second = static_cast<BAM::BaseEnum>( _pickIndexWithSecondHighestMetric(_alleleCounts, majorityBase.get()) );
		_altAlleles.emplace_back(second);
	} else {
		_altAlleles.emplace_back(majorityBase);
		_calledGenotype = "1";
	}

	return true;
};

bool TCallerConsensify::_callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(site.depth() > _downsampleDepth){
		throw std::runtime_error("bool TCallerConsensify::_callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods): depth > _downsampleDepth!");
	}

	//get per allele counts
	_countAlleles(site);

	//call majority
	BAM::Base majorityBase = static_cast<BAM::BaseEnum>( _pickIndexWithHighestMetric(_alleleCounts) );

	//check if we have sufficient depth to call
	if(_alleleCounts[majorityBase.get()] < _minMajorityDepth){
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
TCallerAllelePresence::TCallerAllelePresence(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TCaller(Parameters, Logfile, RandomGenerator){
	//caller settings
	_callerName = "Allele Presence Caller";
	_filenameExtention = "_allelePresence.vcf";
	Logfile->list("Will use the " + _callerName + ".");
	_usesPrior = true;

	//parse VCF fields
	_setAcceptableFields(&_VCFGenotypeFields, "GT,DP,AD,GQ,AP");
	initializeOutput(Parameters, Logfile);

	//initialize allele counts
	MAP = BAM::N;
};

void TCallerAllelePresence::_fillPosteriors(TGenotypeLikelihoods & genotypeLikelihoods){
	//calculate posterior probabilities
	posterior.fill(genotypeLikelihoods, *_genotypePrior);

	//sum for each base
	allelePostProb[0] = posterior[BAM::AA].get() + posterior[BAM::AC].get() + posterior[BAM::AG].get() + posterior[BAM::AT].get();
	allelePostProb[1] = posterior[BAM::AC].get() + posterior[BAM::CC].get() + posterior[BAM::CG].get() + posterior[BAM::CT].get();
	allelePostProb[2] = posterior[BAM::AG].get() + posterior[BAM::CG].get() + posterior[BAM::GG].get() + posterior[BAM::GT].get();
	allelePostProb[3] = posterior[BAM::AT].get() + posterior[BAM::CT].get() + posterior[BAM::GT].get() + posterior[BAM::TT].get();
};

bool TCallerAllelePresence::_callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(!_priorSet) throw "Can not call AllelePresence genotypes: prior has not been set!";

	//fill posteriors for each allele
	_fillPosteriors(genotypeLikelihoods);

	//find MAP
	MAP = static_cast<BAM::BaseEnum>( _pickIndexWithHighestMetric(allelePostProb) );

	//decide on alt
	if(MAP == referenceBase){
		_calledGenotype = "0";

		//find second most common as alternative allele
		BAM::Base second = static_cast<BAM::BaseEnum>( _pickIndexWithSecondHighestMetric(allelePostProb, MAP.get()) );
		_altAlleles.push_back(second);
	} else {
		_altAlleles.push_back(MAP);
		_calledGenotype = "1";
	}

	return true;
};

bool TCallerAllelePresence::_callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(!_priorSet) throw "Can not call AllelePresence genotypes: prior has not been set!";

	//fill posteriors for each allele
	_fillPosteriors(genotypeLikelihoods);

	//find MAP
	uint8_t highest;
	if(allelePostProb[referenceBase.get()] > allelePostProb[_altAlleles[0].get()]){
		highest = 0;
	} else if(allelePostProb[referenceBase.get()] > allelePostProb[_altAlleles[0].get()]){
		highest = 1;
	} else {
		highest = _randomGenerator->sample(2);
	}

	if(!_allowKnownAllelesCallsDifferentFromBestCall){
		if(allelePostProb[highest] < allelePostProb.max()){
			return false;
		}
	}

	//decide on genotype (index 0 is ref base)
	if(highest == 0){
		_calledGenotype = "0";
		MAP = referenceBase;
	} else {
		_calledGenotype = "1";
		MAP = _altAlleles[0];
	}

	return true;
};

std::string TCallerAllelePresence::_getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	return (std::string) BAM::PhredIntErrorRate(allelePostProb[MAP].inverse());
};

std::string TCallerAllelePresence::_getVCFGenotypeString_AP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	std::string ret = (std::string) BAM::PhredIntErrorRate(posterior[0]);
	ret += ',' + (std::string) BAM::PhredIntErrorRate(posterior[1]);
	ret += ',' + (std::string) BAM::PhredIntErrorRate(posterior[2]);
	ret += ',' + (std::string) BAM::PhredIntErrorRate(posterior[3]);
	return ret;
};

/////////////////////////////////////////////////////////
// TCallerDiploid
// common function between MLE, Bayes and gvcf callers
/////////////////////////////////////////////////////////
TCallerDiploid::TCallerDiploid(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TCaller(Parameters, Logfile, RandomGenerator){
	imbalanceCalculated = false;
	_missingGenotype = "./.";

	_setAcceptableFields(&_VCFGenotypeFields, "AB,AI");
};

void TCallerDiploid::_clearAfterCall(){
	TCaller::_clearAfterCall();
	imbalanceCalculated = false;
};

void TCallerDiploid::callGenotypeFromMetric(TGenotypeData & metric){
	genotypeAtMax = static_cast<BAM::GenotypeEnum>( _pickIndexWithHighestMetric(metric) );
	genotypeAtSecond = static_cast<BAM::GenotypeEnum>( _pickIndexWithSecondHighestMetric(metric, genotypeAtMax.get()) );

	//decide on alternative alleles
	if(genotypeAtMax.firstAllele() == referenceBase){
		if(genotypeAtMax.secondAllele() == referenceBase){
			_calledGenotype = "0/0";
			//MLE is homozygous reference -> find second best allele
			if(genotypeAtSecond.firstAllele() == referenceBase)
				_altAlleles.push_back(genotypeAtSecond.secondAllele());
			else if(genotypeAtSecond.secondAllele() == referenceBase)
				_altAlleles.push_back(genotypeAtSecond.firstAllele());
			else {
				//none is ref, pick at random
				_altAlleles.push_back( genotypeAtSecond.randomAllele(*_randomGenerator) );
			}
		} else {
			_altAlleles.push_back(genotypeAtMax.secondAllele());
			_calledGenotype = "0/1";
		}
	} else {
		if(genotypeAtMax.secondAllele() == referenceBase){
			_altAlleles.push_back(genotypeAtMax.firstAllele());
			_calledGenotype = "0/1";
		} else {
			if(genotypeAtMax.isHomozygous()){
				_altAlleles.push_back(genotypeAtMax.firstAllele());
				_calledGenotype = "1/1";

				//find second best allele, but give preference to reference in case likelihoods are equal
				if(_allowTriallelicSites){
					//int hetRef = genoMap.getGenotype(referenceBase, genoMap.genotypeToBase[indexOfMax][0]);
					//int homRef = genoMap.getGenotype(referenceBase, referenceBase);

					//only use second alternative allele in case het genotype with reference is less likely

					if(referenceBase == BAM::N ||
					   (metric[ Genotype(referenceBase, genotypeAtMax.firstAllele())] < metric[ genotypeAtSecond]
						&& metric[ Genotype(referenceBase, referenceBase) ] < metric[ genotypeAtSecond ])){
						if(genotypeAtSecond.firstAllele() == referenceBase || genotypeAtSecond.firstAllele() == _altAlleles[0])
							_altAlleles.push_back(genotypeAtSecond.secondAllele());
						else if(genotypeAtSecond.secondAllele() == referenceBase || genotypeAtSecond.secondAllele() == _altAlleles[0])
							_altAlleles.push_back(genotypeAtSecond.firstAllele());
						else {
							//none is ref, pick at random
							_altAlleles.push_back( genotypeAtSecond.randomAllele(*_randomGenerator) );
						}
					}
				}
			} else {
				if(_allowTriallelicSites){
					//allow triallelic sites
					_altAlleles.push_back(genotypeAtMax.firstAllele());
					_altAlleles.push_back(genotypeAtMax.secondAllele());
					_calledGenotype = "1/2";
				} else {
					//decide on which of the two alternative alleles to pick -> check second highest
					if(genotypeAtSecond.isHomozygous()){
						if(genotypeAtSecond.firstAllele() == genotypeAtMax.firstAllele())
							_altAlleles.push_back(genotypeAtMax.firstAllele());
						else if(genotypeAtSecond.firstAllele() == genotypeAtMax.secondAllele())
							_altAlleles.push_back(genotypeAtMax.secondAllele());
						else {
							//neither alt 0 nor 1, pick at random
							_altAlleles.push_back(genotypeAtMax.randomAllele(*_randomGenerator));
						}
					} else {
						if(genotypeAtSecond.firstAllele() == referenceBase && (genotypeAtSecond.secondAllele() == genotypeAtMax.firstAllele() || genotypeAtSecond.secondAllele() == genotypeAtMax.firstAllele())){
							_altAlleles.push_back(genotypeAtSecond.secondAllele());
						} else if(genotypeAtSecond.secondAllele() == referenceBase && (genotypeAtSecond.firstAllele() == genotypeAtMax.firstAllele() || genotypeAtSecond.firstAllele() == genotypeAtMax.secondAllele())){
							_altAlleles.push_back(genotypeAtSecond.firstAllele());
						} else {
							//neither alt 0 nor 1, pick at random
							_altAlleles.push_back(genotypeAtMax.randomAllele(*_randomGenerator));
						}
					}

					//now call genotype from these alleles
					callGenotypeFromMetricKnownAlleles(metric);
				}
			}
		}
	}
};

void TCallerDiploid::callGenotypeFromMetricKnownAlleles(const TGenotypeData & metric){
	//get genotypes
	BAM::BiallelicGenotypes geno(referenceBase, _altAlleles[0]);

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

	_calledGenotype = vec[_randomGenerator->sample(vec.size())];
};

bool TCallerDiploid::callGenotypeFromMetricKnownAllelesUpdateIndex(const TGenotypeData & metric){
	//initialize
	static std::vector<std::string> gt = {"0/0", "0/1", "1/1"};

	//get genotypes
	BAM::BiallelicGenotypes geno(referenceBase, _altAlleles[0]);

	//subset metric to considered genotypes
	std::array<double, 3> metricKnownAlleles;
	metricKnownAlleles[0] = metric[geno.homoFirst()];
	metricKnownAlleles[1] = metric[geno.het()];
	metricKnownAlleles[2] = metric[geno.homoSecond()];

	uint8_t best = _pickIndexWithHighestMetric(metricKnownAlleles);
	uint8_t secondBest = _pickIndexWithSecondHighestMetric(metricKnownAlleles, best);

	genotypeAtMax = static_cast<BAM::GenotypeEnum>(geno[best] );
	genotypeAtSecond = static_cast<BAM::GenotypeEnum>( geno[secondBest] );
	_calledGenotype = gt[best];

	if(!_allowKnownAllelesCallsDifferentFromBestCall){
		//check if call matches metric of best call
		if(metric[genotypeAtMax] < metric.max()){
			return false;
		}
	}
	return true;
};

std::string TCallerDiploid::getPerGenotypeMetricString(TGenotypeData & metric){
	//if you have alleles R, A, B, C then the order of the PL is: RR, RA, AA | RB, AB, BB | RC, AC, BC, CC
	//plot missing value (.) for all metrics involving the reference if the reference is N
	std::string ret;
	//first for reference base
	if(referenceBase == BAM::N)
		ret = ".";
	else
		ret = toString(metric[Genotype(referenceBase, referenceBase)]);

	//now for alternative alleles
	if(_altAlleles.size() > 0){
		if(referenceBase == BAM::N)
			ret += ",.";
		else
			ret += ',' + toString(metric[Genotype(referenceBase, _altAlleles[0])]);
		ret += ',' + toString(metric[Genotype(_altAlleles[0], _altAlleles[0])]);

		if(_altAlleles.size() > 1){
			if(referenceBase == BAM::N)
				ret += ",.";
			else
				ret += ',' + toString(metric[Genotype(referenceBase, _altAlleles[1])]);
			ret += ',' + toString(metric[Genotype(_altAlleles[0], _altAlleles[1])]);
			ret += ',' + toString(metric[Genotype(_altAlleles[1], _altAlleles[1])]);
		}

		if(_altAlleles.size() > 2){
			if(referenceBase == BAM::N)
				ret += ",.";
			else
				ret += ',' + toString(metric[Genotype(referenceBase, _altAlleles[2])]);
			ret += ',' + toString(metric[Genotype(_altAlleles[0], _altAlleles[2])]);
			ret += ',' + toString(metric[Genotype(_altAlleles[1], _altAlleles[2])]);
			ret += ',' + toString(metric[Genotype(_altAlleles[2], _altAlleles[2])]);
		}
	}
	return ret;
};

void TCallerDiploid::calculateImbalance(const TSite & site){
	if(!imbalanceCalculated){
		if(!_altAlleles.empty()){
			_countAlleles(site);

			if(_altAlleles.size() == 1){
				double sum = (_alleleCounts[referenceBase] + _alleleCounts[_altAlleles[0]]);
				if(referenceBase == BAM::N || sum == 0){
					AB = '.'; AI = '.';
				} else {
					AB = toString(_alleleCounts[referenceBase] / sum);
					AI = toString(_binomP.binomPValue(_alleleCounts[referenceBase], _alleleCounts[_altAlleles[0]]));
				}
			} else {
				if(genotypeAtMax.isHeterozygous()){
					double sum = (double) _alleleCounts[genotypeAtMax.firstAllele()] + _alleleCounts[genotypeAtMax.secondAllele()];
					if(sum == 0){
						AB = '.'; AI = '.';
					} else {
						AB = toString(_alleleCounts[genotypeAtMax.firstAllele()] / sum);
						AI = toString(_binomP.binomPValue(_alleleCounts[genotypeAtMax.firstAllele()], _alleleCounts[genotypeAtMax.secondAllele()]));
					}
				} else { // is homo -> do it against the second alternative allele
					double sum = (double) _alleleCounts[_altAlleles[0]] + _alleleCounts[_altAlleles[1]];
					if(sum == 0){
						AB = '.'; AI = '.';
					} else {
						AB = toString(_alleleCounts[_altAlleles[0]] / sum);
						AI = toString(_binomP.binomPValue(_alleleCounts[_altAlleles[0]], _alleleCounts[_altAlleles[1]]));
					}
				}
			}
		}
		imbalanceCalculated = true;
	}
};

std::string TCallerDiploid::_getVCFGenotypeString_AB(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(_altAlleles.empty()) return ".";

	calculateImbalance(site);
	return AB;
};

std::string TCallerDiploid::_getVCFGenotypeString_AI(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(_altAlleles.empty()) return ".";

	calculateImbalance(site);
	return AI;
};

/////////////////////////////////////////////////////////
// TCallerMLE
/////////////////////////////////////////////////////////
TCallerMLE::TCallerMLE(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TCallerDiploid(Parameters, Logfile, RandomGenerator){
	//caller settings
	_callerName = "MLE Caller";
	_filenameExtention = "_maximumLikelihood.vcf";
	Logfile->list("Will use the " + _callerName + ".");

	//parse VCF fields
	_setAcceptableFields(&_VCFGenotypeFields, "AD,GQ,GL,PL");
	printGenotypeFields("GT,DP,AD,GQ,PL"); //set default tags to print
	initializeOutput(Parameters, Logfile);
};

bool TCallerMLE::_callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	callGenotypeFromMetric(genotypeLikelihoods);
	return true;
};

bool TCallerMLE::_callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	return callGenotypeFromMetricKnownAllelesUpdateIndex(genotypeLikelihoods);
};

std::string TCallerMLE::_getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	Probability error(genotypeLikelihoods[genotypeAtSecond].get() / genotypeLikelihoods[genotypeAtMax].get());
	BAM::PhredIntErrorRate phred(error);
	return toString(phred);
};

std::string TCallerMLE::_getVCFGenotypeString_GL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//normalize
	for(BAM::Genotype g = BAM::Genotype::min(); g < BAM::Genotype::max(); ++g){
		tmpGenoData[g] = log10(genotypeLikelihoods[g].get() / genotypeLikelihoods[genotypeAtMax].get());
	}

	//get string
	return getPerGenotypeMetricString(tmpGenoData);
};

std::string TCallerMLE::_getVCFGenotypeString_PL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//normalize
	double phredMax = _qualMap.errorToPhred(genotypeLikelihoods[indexOfMax]);
	for(uint8_t g=0; g<10; ++g)
		tmpGenoData[g] = (int) round(_qualMap.errorToPhred(genotypeLikelihoods[g]) - phredMax);

	//get string
	return getPerGenotypeMetricString(tmpGenoData);
};

//------------------------------------------------------
// TCallerBayes
//------------------------------------------------------
TCallerBayes::TCallerBayes(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TCallerDiploid(Parameters, Logfile, RandomGenerator){
	//caller settings
	_callerName = "Bayesian Caller";
	_filenameExtention = "_maximumAPosteriori.vcf";
	_usesPrior = true;
	Logfile->list("Will use the " + _callerName + ".");

	//parse VCF fields
	_setAcceptableFields(&_VCFGenotypeFields, "AD,GQ,GP");
	printGenotypeFields("GT,DP,AD,GQ,GP");  //set default tags to print
	initializeOutput(Parameters, Logfile);
};


bool TCallerBayes::_callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(!_priorSet) throw "Can not call Bayesian genotypes: prior has not been set!";

	//calculate posterior probabilities
	posterior.fill(genotypeLikelihoods, *_genotypePrior);

	//call
	callGenotypeFromMetric(posterior);
	return true;
};

bool TCallerBayes::_callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(!_priorSet) throw "Can not call Bayesian genotypes: prior has not been set!";

	//calculate posterior probabilities
	posterior.fill(genotypeLikelihoods, *_genotypePrior);

	//call
	return callGenotypeFromMetricKnownAllelesUpdateIndex(posterior);
};

std::string TCallerBayes::_getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	return toString(_qualMap.errorToPhredInt(1.0 - posterior[indexOfMax]));
};

std::string TCallerBayes::_getVCFGenotypeString_GP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//phred
	for(uint8_t g=0; g<10; ++g)
		tmpGenoData[g] = (int) _qualMap.errorToPhredInt(posterior[g]);

	//get string
	return getPerGenotypeMetricString(tmpGenoData);
};

//------------------------------------------------------
// TCall
// the class to perform calls based on windows
//------------------------------------------------------
TCall::TCall(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_windows(Parameters, Logfile, RandomGenerator){
	//initialize caller
	_logfile->startIndent("Initializing caller:");
	std::string method = Parameters.getParameterWithDefault<std::string>("method", "MLE");
	if(method == "randomBase"){
		_caller = new TCallerRandomBase(Parameters, _logfile, _randomGenerator);
	} else if(method == "majorityBase"){
		_caller = new TCallerMajorityBase(Parameters, _logfile, _randomGenerator);
	} else if(method == "consensify"){
		_caller = new TCallerConsensify(_downsampleDepth, Parameters, _logfile, _randomGenerator);
	} else if(method == "allelePresence"){
		_caller = new TCallerAllelePresence(Parameters, _logfile, _randomGenerator);
	} else if(method == "MLE"){
		_caller = new TCallerMLE(Parameters, _logfile, _randomGenerator);
	} else if(method == "Bayesian"){
		_caller = new TCallerBayes(Parameters, _logfile, _randomGenerator);
	} else if(method == "gVCF"){
		throw "GVCF NOT YET IMPLEMENTED!";
		_caller->printSitesWithNoData();
	} else throw "Unknown calling method '" + method + "'! Use randomBase, allelePresence, MLE, Bayesian or gVCF.";
	_logfile->list("Will use the " + _caller->name() + ".");

	//prior setting
	if(_caller->usesPrior()){
		_initializeGenotypePrior(Parameters);
	} else {
		_prior = new TGenotypePrior();
	}
	_caller->setPrior(_prior->getPointerToPrior());

	//read output settings
	_caller->initializeOutput(Parameters, _logfile);

	//open output file
	std::string sampleName = Parameters.getParameterWithDefault<std::string>("indName", _outputName);
	_logfile->list("Will use sample name '" + sampleName + "'. (parameter 'sampleName')");
	_caller->openVCF(_outputName + "_calls", sampleName, _logfile);

	//limit to sites with known alleles?
	if(Parameters.parameterExists("alleles")){
		_logfile->startIndent("Will limit calls to sites with known alleles (parameter 'alleles'):");
		_openSiteSubset("alleles");
		_logfile->endIndent();
	} else {
		_logfile->list("Will call without prior knowledge on alleles. (use 'alleles' to provide known alleles)");
		//make sure FASTA is open unless alleles are provided
		_openReference(true);
	}
	_logfile->endIndent();
};

TCall::~TCall(){
	delete _caller;
	delete _prior;
};

void TCall::_initializeGenotypePrior(TParameters & Parameters){
	_logfile->startIndent("Initializing genotype prior:");
	//read prior from parameters
	std::string priorMethod = Parameters.getParameterWithDefault<std::string>("prior", "theta");
	if(priorMethod == "unif"){
		_prior = new TGenotypePriorUniform();
		_logfile->list("Will use a uniform prior with equal weights for all genotypes.");
	} else if(priorMethod == "theta"){
		if(Parameters.parameterExists("fixedTheta")){
			double theta = Parameters.getParameter<double>("fixedTheta");
			_logfile->list("Will use a fixed theta = " + toString(theta));
			bool equalBaseFreq = Parameters.parameterExists("equalBaseFreq");
			if(equalBaseFreq)
				_logfile->list("Will use equal base frequencies.");
			else
				_logfile->list("Will estimate base frequencies individually for each window.");
			_prior = new TGenotypePriorFixedTheta(theta, equalBaseFreq, _logfile, _randomGenerator);
		} else {
			_logfile->list("Will use a prior based on theta and base frequencies estimated individually for each window.");
			std::string thetaOuputName = _outputName + "_theta_estimates.txt.gz";
			if(Parameters.parameterExists("defaultTheta")){
				double defaultTheta = Parameters.getParameter<double>("defaultTheta");
				_logfile->list("Will use a default theta of ", defaultTheta, " for windows with limited data.");
				_prior = new TGenotypePriorTheta(Parameters, thetaOuputName, defaultTheta, _logfile, _randomGenerator);
			} else
				_prior = new TGenotypePriorTheta(Parameters, thetaOuputName, _logfile, _randomGenerator);
		}
	} else throw "Unknown prior type '" + priorMethod + "'!";
	_logfile->endIndent();
};

void TCall::_call(){
	uint32_t pos = 0;
	for(auto& s : _window){
		_genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(s, _genoLik);
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
		std::set<TSiteSubsetSite> thesePositions = _subset->getPositionInWindow(_window);
		for(auto& it : thesePositions){
			//calculate genotype likelihoods
			uint32_t internalPos = it - _window.from();
			TSite& site = _window[internalPos];
			site.setRefBase(it.ref());
			_genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(site, _genoLik);
			_caller->call(_window.chrName(), _window.positionOnChr(internalPos), site, _genoLik, it.ref(), it.alt());
		}
	}
};


void TCall::_handleWindow(){
	if(_window.passedFilters() || _caller->printSitesWithNoData()){
		//update genotype prior
		GenotypeLikelihoods::TBaseData freq;
		_window.estimateBaseFrequencies(freq);
		_prior->update(_window, _logfile, _genotypeLikelihoodCalculator);

		//call
		_logfile->listFlushTime("Calling genotypes ...");
		if(_subset){
			_callKnwonAlleles();
		} else {
			_call();
		}
		_logfile->doneTime();
	}
};

void TCall::call(){
	_traverseBAMWindows();
};

}; // end namespace
