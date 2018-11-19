/*
 * TCaller.cpp
 *
 *  Created on: Nov 17, 2018
 *      Author: phaentu
 */


#include "TCaller.h"

/////////////////////////////////////////////////////////
// TCaller
/////////////////////////////////////////////////////////
TCaller::TCaller(TRandomGenerator* RandomGenerator){
	callerName = "default caller";
	filenameExtention = ".vcf.gz";
	randomGenerator = RandomGenerator;

	//output choices
	_printSitesWithNoData = false;

	//vcf file
	vcfOpen = false;

	//set acceptable tags
	setAcceptableFields(&VCFInfoFields, "DP");
	setAcceptableFields(&VCFGenotypeFields, "GT,DP");

	//set default tags to print
	printInfoFields("DP");
	printGenotypeFields("GT,DP");

	//tmp variables
	numGenotypes = 10;
};

TCaller::~TCaller(){
	closeVCF();
};

//-------------------------------------------------------------------------------------------
// Output settings
//-------------------------------------------------------------------------------------------
void TCaller::setAcceptableFields(TVCFFieldVector* fields, std::string tags){
	std::vector<std::string> vec;
	fillVectorFromStringAnySkipEmpty(tags, vec, ",");
	for(std::vector<std::string>::iterator it = vec.begin(); it != vec.end(); ++it)
		fields->acceptField(*it);
};

void TCaller::printField(TVCFFieldVector* fields, std::string tag){
	if(!fields->useField(tag))
		throw "VCF " +  fields->type() + " field '" + tag + "' can not be printed by the " + callerName +"!";
};

void TCaller::printInfoFields(std::vector<std::string> & tags){
	VCFInfoFields.clearUsed();
	for(std::vector<std::string>::iterator it = tags.begin(); it != tags.end(); ++it){
		printField(&VCFInfoFields, *it);
	}
	fillInfoFieldFunctionPointers();
};

void TCaller::printInfoFields(std::string tags){
	std::vector<std::string> vec;
	fillVectorFromStringAnySkipEmpty(tags, vec, ",");
	printInfoFields(vec);
};

void TCaller::printGenotypeFields(std::vector<std::string> & tags){
	VCFGenotypeFields.clearUsed();
	for(std::vector<std::string>::iterator it = tags.begin(); it != tags.end(); ++it)
		printField(&VCFGenotypeFields, *it);
	fillGenotypeFieldFunctionPointers();
};

void TCaller::printGenotypeFields(std::string tags){
	std::vector<std::string> vec;
	fillVectorFromStringAnySkipEmpty(tags, vec, ",");
	printGenotypeFields(vec);
};

void TCaller::reportSettings(TLog* logfile){
	//report caller name
	logfile->list("Will use the " + callerName + ".");

	//report VCF fields
	VCFInfoFields.reportUsedFields(logfile);
	VCFGenotypeFields.reportUsedFields(logfile);

	//report whether all sites are printed
	if(_printSitesWithNoData)
		logfile->list("Will print all sites, also those without data");
}

//-------------------------------------------------------------------------------------------
// open / close VCF file, print header
//-------------------------------------------------------------------------------------------
void TCaller::openVCF(const std::string FilenameTag, const std::string sampleName){
	filename = FilenameTag  + filenameExtention + ".gz";
	vcf.open(filename.c_str());
	if(!vcf) throw "Failed to open VCF file '" + filename + "' for writing!";
	vcfOpen = true;

	//write header
	writeVCFHeader(sampleName);
};

void TCaller::closeVCF(){
	if(vcfOpen){
		vcf.close();
		vcfOpen = false;
	}
};

void TCaller::writeVCFHeader(const std::string & sampleName){
	//write header
	vcf << "##fileformat=VCFv4.2\n";
	vcf << "##source=atlas\n";

	//write INFO and GENOTYPE fields
	VCFInfoFields.writeVCFHeader(vcf);
	VCFGenotypeFields.writeVCFHeader(vcf);

	//write column header
	vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleName << "\n";
};



void TCaller::callGenotype(TSite & site){


	calledGenotype = '0';
	std::vector<int> altAlleles;



};

//-------------------------------------------------------------------------------------------
// Info fields
//-------------------------------------------------------------------------------------------
void TCaller::fillInfoFieldFunctionPointers(){
	//clear current vector
	VCFInfoFunctionsVec.clear();

	//get used tags
	std::vector<std::string> tagVec;
	VCFInfoFields.fillVectorWithTagsOfUsedFields(tagVec);

	//add functions to info field vector
	for(std::vector<std::string>::iterator it = tagVec.begin(); it != tagVec.end(); ++it){
		if(*it == "DP")
			VCFInfoFunctionsVec.push_back( &TCaller::getVCFInfoString_DP );
		else throw "No function defined for VCF " + VCFInfoFields.type() + " field '" + *it + "'! @Programmer: add function to TTCaller::fillInfoFieldFunctionPointers()!";
	}

};

std::string TCaller::getVCFInfoString_DP(const TSite & site){
	std::cout << "==============> IN FUNCTION TCaller::getVCGInfoString_DP(const TSite & site)" << std::endl;
	return "DP=" + toString(site.bases.size());
};

//-------------------------------------------------------------------------------------------
// genotype fields
//-------------------------------------------------------------------------------------------
void TCaller::fillGenotypeFieldFunctionPointers(){
	//clear current vector
	VCFGenotypeFunctionsVec.clear();
	genotypeFormatString.clear();

	//get used tags
	std::vector<std::string> tagVec;
	VCFGenotypeFields.fillVectorWithTagsOfUsedFields(tagVec);

	//add functions to genotype field vector
	for(std::vector<std::string>::iterator it = tagVec.begin(); it != tagVec.end(); ++it){
		if(*it == "GT")
			VCFGenotypeFunctionsVec.push_back( &TCaller::getVCFGenotypeString_GT );
		else if(*it == "DP")
			VCFGenotypeFunctionsVec.push_back( &TCaller::getVCFGenotypeString_DP );
		else if(*it == "GQ")
			VCFGenotypeFunctionsVec.push_back( &TCaller::getVCFGenotypeString_GQ );
		else throw "No function defined for VCF " + VCFGenotypeFields.type() + " field '" + *it + "'! @Programmer: add function to TTCaller::fillGenotypeFieldFunctionPointers()!";

		//add to format string
		if(genotypeFormatString.length() > 0) genotypeFormatString += ':';
		genotypeFormatString += *it;
	}
};

std::string TCaller::getVCFGenotypeString_GT(const TSite & site){
	return calledGenotype;
};

std::string TCaller::getVCFGenotypeString_DP(const TSite & site){
	return toString(site.bases.size());
};


//-------------------------------------------------------------------------------------------
// writing VCF
//-------------------------------------------------------------------------------------------
std::string TCaller::composeVCFString(std::vector<std::string (TCaller::*)(const TSite & site)> & vec, const TSite & site){
	//no info fields?
	if(vec.empty()) return ".";

	//add first info
	std::vector<std::string (TCaller::*)(const TSite & site)>::iterator it = vec.begin();
	std::string info = (this->*(*it))(site);
	++it;

	//loop over res
	for(; it != vec.end(); ++it)
		info += ';' + (this->*(*it))(site);

	return info;
};

void TCaller::writeCallToVCF(const std::string & chr, const long pos, TSite & site){
	//write chr, position and (no) variant ID
	vcf << chr << '\t' << pos << "\t.\t";

	//write reference and alternative alleles
	vcf << site.referenceBase << "\t";
	if(altAlleles.size() == 0){
		vcf << '.';
	} else {
		vcf << altAlleles[0];
		for(size_t i=1; i<altAlleles.size(); ++i)
			vcf << ',' << altAlleles[i];
	}

	//write (no) variant quality and (no) filter
	vcf << "\t.\t.";

	//write info fields
	vcf << '\t' << composeVCFString(VCFInfoFunctionsVec, site);

	//write genotype fields
	vcf << '\t' << genotypeFormatString << '\t' << composeVCFString(VCFGenotypeFunctionsVec, site);

	//end with new line
	vcf << '\n';

	//clean up storage
	altAlleles.clear();
};

/*
void TSite::calculateNormalizedGenotypeLikelihoodsAndQuality(TRandomGenerator & randomGenerator, double* emissionProbabilitiesPhredScaled,  double & quality, double & maxGenotypeProb, int & MLGenotype){
	//calculate phred-scaled likelihoods and find max
	maxGenotypeProb = 100000.0;
	quality = 100000.0;
	std::vector<int> MLEs;
	std::vector<int>::iterator it;
	for(int i=0; i<numGenotypes; ++i){
		emissionProbabilitiesPhredScaled[i] = makePhredByRef(emissionProbabilities[i]);
		if(emissionProbabilitiesPhredScaled[i] < maxGenotypeProb){
			MLGenotype = i;
			quality = maxGenotypeProb;
			maxGenotypeProb = emissionProbabilitiesPhredScaled[i];
			MLEs.clear();
			MLEs.push_back(i);
		} else if(emissionProbabilitiesPhredScaled[i] == maxGenotypeProb){
			MLEs.push_back(i);
			quality = emissionProbabilitiesPhredScaled[i];
		} else if(emissionProbabilitiesPhredScaled[i] < quality){
			quality = emissionProbabilitiesPhredScaled[i];
		}
	}
	//select best allele at random if there are multiple options
	MLGenotype = MLEs[randomGenerator.pickOne(MLEs.size())];
	quality = quality - maxGenotypeProb;
}
*/


void TCaller::call(const std::string & chr, const long pos, TSite & site){
	//check if there is data
	if(site.hasData){
		//call
		callGenotype(site);

		//check if we write
		writeCallToVCF(chr, pos, site);

	} else {
		if(_printSitesWithNoData)
			vcf << "\t.\t" << site.referenceBase << "\t.\t.\t.\t.\tGT:DP:GQ\t./.:0:0"; //TO: check for Bayesian case?
	}
};


/////////////////////////////////////////////////////////
// TCallerRandomBase
/////////////////////////////////////////////////////////
TCallerRandomBase::TCallerRandomBase(TRandomGenerator* RandomGenerator):TCaller(RandomGenerator){
	//caller settings
	callerName = "Random Base Caller";
	filenameExtention = "_randomBase.vcf";
	defaultInfoFields = "DP";
	defaultGenotypeFields = "GT,DP";
};

void TCallerRandomBase::callGenotype(TSite & site){
	//randomly pick a base
	char allele = site.bases[randomGenerator->pickOne(site.bases.size())]->getBase();

	//decide on alt
	if(allele == site.referenceBase){
		calledGenotype = "0";
	} else {
		altAlleles.push_back(allele);
		calledGenotype = "1";
	}
};


/////////////////////////////////////////////////////////
// TCallerDiploid
/////////////////////////////////////////////////////////
TCallerDiploid::TCallerDiploid(TRandomGenerator* RandomGenerator):TCaller(RandomGenerator){
	peformImbalanceTest = false;
};

//------------------------------------------------------
// TCallerMLE
//------------------------------------------------------
TCallerMLE::TCallerMLE(TRandomGenerator* RandomGenerator):TCallerDiploid(RandomGenerator){

};

void TCallerMLE::fillCallingMetric(TSite & site){
	//fill with likelihoods, normalized by maximum
	//find maximum
	double maxGenotypeProb = site.emissionProbabilities[0];
	for(int i=1; i<numGenotypes; ++i){
		if(site.emissionProbabilities[i] > maxGenotypeProb)
			maxGenotypeProb = site.emissionProbabilities[i];
	}

	//normalize
	for(int i=0; i<numGenotypes; ++i){
		perGenotypeMetric[i] = site.emissionProbabilities[i] / maxGenotypeProb;
	}
};



//------------------------------------------------------
// TCallerBayes
//------------------------------------------------------
TCallerBayes::TCallerBayes(TRandomGenerator* RandomGenerator):TCaller(RandomGenerator){

};

void TCallerBayes::setPrior(double* GenotypePriorProbabilities){
	for(int i=0; i<numGenotypes; ++i)
		genotypePriorProbabilities[i] = GenotypePriorProbabilities[i];
};

void TCallerBayes::fillCallingMetric(TSite & site){
	double tot = 0.0;

	for(int i=0; i<numGenotypes; ++i){
		perGenotypeMetric[i] = site.emissionProbabilities[i] * genotypePriorProbabilities[i];
		tot += perGenotypeMetric[i];
	}

	//normalize
	for(int i=0; i<numGenotypes; ++i)
		perGenotypeMetric[i] /= tot;
};

