/*
 * TCaller.cpp
 *
 *  Created on: Nov 17, 2018
 *      Author: phaentu
 */


#include "TCaller.h"

using namespace GenotypeLikelihoods;

/////////////////////////////////////////////////////////
// TCaller
/////////////////////////////////////////////////////////
TCaller::TCaller(TRandomGenerator* RandomGenerator){
	callerName = "default caller";
	filenameExtention = ".vcf.gz";
	randomGenerator = RandomGenerator;

	//output choices
	_printSitesWithNoData = false;
	_printAltIfHomoRef = true;
	_allowTriallelicSites = true;
	missingGenotype = ".";

	//vcf file
	vcfOpen = false;

	//set acceptable tags
	setAcceptableFields(&VCFInfoFields, "DP");
	setAcceptableFields(&VCFGenotypeFields, "GT,DP,AD");
	_usesPrior = false;
	priorSet = false;

	//set default tags to print
	printInfoFields("");
	printGenotypeFields("GT,DP");

	//tmp variables
	allelesCounted = false;
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

std::string TCaller::getVCFInfoString_DP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
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
		else if(*it == "AD")
			VCFGenotypeFunctionsVec.push_back( &TCaller::getVCFGenotypeString_AD );
		else if(*it == "AP")
			VCFGenotypeFunctionsVec.push_back( &TCaller::getVCFGenotypeString_AP );
		else if(*it == "GL")
			VCFGenotypeFunctionsVec.push_back( &TCaller::getVCFGenotypeString_GL );
		else if(*it == "PL")
			VCFGenotypeFunctionsVec.push_back( &TCaller::getVCFGenotypeString_PL );
		else if(*it == "GP")
			VCFGenotypeFunctionsVec.push_back( &TCaller::getVCFGenotypeString_GP );
		else if(*it == "AB")
			VCFGenotypeFunctionsVec.push_back( &TCaller::getVCFGenotypeString_AB );
		else if(*it == "AI")
			VCFGenotypeFunctionsVec.push_back( &TCaller::getVCFGenotypeString_AI );
		else throw "No function defined for VCF " + VCFGenotypeFields.type() + " field '" + *it + "'! @Programmer: add function to TTCaller::fillGenotypeFieldFunctionPointers()!";

		//add to format string
		if(genotypeFormatString.length() > 0) genotypeFormatString += ':';
		genotypeFormatString += *it;
	}
};

std::string TCaller::getVCFGenotypeString_GT(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	return calledGenotype;
};

std::string TCaller::getVCFGenotypeString_DP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	return toString(site.bases.size());
};

std::string TCaller::getVCFGenotypeString_AD(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	countAlleles(site);
	std::string ret;
	if(referenceBase == N) ret = "0";
	else ret = toString(alleleCounts[referenceBase]);

	for(std::vector<int>::iterator it = altAlleles.begin(); it != altAlleles.end(); ++it)
		ret += ',' + toString(alleleCounts[*it]);
	return ret;
};


//-------------------------------------------------------------------------------------------
// writing VCF
//-------------------------------------------------------------------------------------------
std::string TCaller::composeVCFString(std::vector<std::string (TCaller::*)(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods)> & vec, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
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

void TCaller::writeAlternativeAllelesToVCF(){
	if(altAlleles.size() == 0){
		vcf << '.';
	} else {
		vcf << genoMap.baseToChar[altAlleles[0]];
		for(size_t i=1; i<altAlleles.size(); ++i)
			vcf << ',' << genoMap.baseToChar[altAlleles[i]];
	}
};

void TCaller::writeCallToVCF(const std::string & chr, const long pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//apply filter on alternative alleles
	if(!_printAltIfHomoRef && (calledGenotype == "0/0" || calledGenotype == "0"))
		altAlleles.clear();

	//write chr, position and (no) variant ID
	vcf << chr << '\t' << pos + 1 << "\t.\t"; //all internal positions are zero-based!

	//write reference and alternative alleles
	vcf << site.referenceBase << "\t";
	writeAlternativeAllelesToVCF();

	//write (no) variant quality and (no) filter
	vcf << "\t.\t.";

	//write info fields
	vcf << '\t' << composeVCFString(VCFInfoFunctionsVec, site, genotypeLikelihoods);

	//write genotype fields
	vcf << '\t' << genotypeFormatString << '\t' << composeVCFString(VCFGenotypeFunctionsVec, site, genotypeLikelihoods);

	//end with new line
	vcf << '\n';

	//clean up storage
	clearAfterCall();
};

void TCaller::writeMissingDataToVCF(const TSite & site){
	if(_printSitesWithNoData)
		vcf << "\t.\t" << site.referenceBase << "\t.\t.\t.\t.\tGT:DP\t" << missingGenotype << ":0";
};

void TCaller::clearAfterCall(){
	altAlleles.clear();
	allelesCounted = false;
};

//-------------------------------------------------------------------------------------------
// calling
//-------------------------------------------------------------------------------------------
void TCaller::countAlleles(const TSite & site){
	if(!allelesCounted){
		site.countAlleles(alleleCounts);
		allelesCounted = true;
	}
};

void TCaller::callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	calledGenotype = "./.";
};

void TCaller::callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	calledGenotype = "./.";
};

void TCaller::call(const std::string & chr, const long pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//set reference base from site
	referenceBase = genoMap.getBase(site.referenceBase);

	//check if there is data
	if(!site.hasData || (referenceBase == N && !_allowTriallelicSites))
		writeMissingDataToVCF(site);
	else {
		//call
		callGenotype(site, genotypeLikelihoods);

		//check if we write
		writeCallToVCF(chr, pos, site, genotypeLikelihoods);
	}
};

void TCaller::call(const std::string & chr, const long pos, const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods, const char & firstAllele, const char & secondAllele){
	//check if there is data
	if(site.hasData){
		//set reference base from site
		referenceBase = genoMap.getBase(site.referenceBase);

		//call
		if(referenceBase == genoMap.getBase(firstAllele))
			altAlleles.push_back(genoMap.getBase(secondAllele));
		else
			altAlleles.push_back(genoMap.getBase(firstAllele));
		callGenotypeKnownAlleles(site, genotypeLikelihoods);

		//check if we write
		writeCallToVCF(chr, pos, site, genotypeLikelihoods);

	} else
		writeMissingDataToVCF(site);
};

template <typename T> int TCaller::pickIndexWithHighestMetric(T* metric, const int size){
	//find maximum
	double maxMetric = 0.0;
	for(int i=0; i<size; ++i){
		if(metric[i] > maxMetric)
			maxMetric = metric[i];
	}

	//get vec of all index at maximum
	std::vector<int> vec;
	for(int i=0; i<size; ++i){
		if(metric[i] == maxMetric)
			vec.push_back(i);
	}

	//return random index among those at max
	return vec[randomGenerator->pickOne(vec.size())];
};

template <typename T> int TCaller::pickIndexWithSecondHighestMetric(T* metric, const int size, const int excludeIndex){
	//find maximum
	double max = 0.0;
	for(int i=0; i<size; ++i){
		if(i != excludeIndex && metric[i] > max)
			max = metric[i];
	}

	//get vec of all index at maximum
	std::vector<int> vec;
	for(int i=0; i<size; ++i){
		if(i != excludeIndex && metric[i] == max)
			vec.push_back(i);
	}

	//return random index among those at max
	return vec[randomGenerator->pickOne(vec.size())];
};


/////////////////////////////////////////////////////////
// TCallerRandomBase
/////////////////////////////////////////////////////////
TCallerRandomBase::TCallerRandomBase(TRandomGenerator* RandomGenerator):TCaller(RandomGenerator){
	//caller settings
	callerName = "Random Base Caller";
	filenameExtention = "_randomBase.vcf";

	//set acceptable tags
	setAcceptableFields(&VCFInfoFields, "DP");
	setAcceptableFields(&VCFGenotypeFields, "GT,DP,AD");
};

void TCallerRandomBase::callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//randomly pick a base
	int allele = site.bases[randomGenerator->pickOne(site.bases.size())]->base;

	//decide on alt
	if(allele == referenceBase){
		calledGenotype = "0";
	} else {
		altAlleles.push_back(allele);
		calledGenotype = "1";
	}
};

void TCallerRandomBase::callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//randomly pick a base among known alleles
	countAlleles(site);
	double probRef = (double) alleleCounts[referenceBase] / (double) (alleleCounts[referenceBase] + alleleCounts[altAlleles[0]]);

	//pick among known alleles
	if(randomGenerator->getRand() < probRef){
		calledGenotype = "0";
	} else {
		calledGenotype = "1";
	}
};

/////////////////////////////////////////////////////////
// TCallerMajorityCall
/////////////////////////////////////////////////////////
TCallerMajorityBase::TCallerMajorityBase(TRandomGenerator* RandomGenerator):TCaller(RandomGenerator){
	//caller settings
	callerName = "Majority Base Caller";
	filenameExtention = "_majorityBase.vcf";
};

void TCallerMajorityBase::callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//get per allele counts
	countAlleles(site);
	int majorityIndex = pickIndexWithHighestMetric(alleleCounts, 4);

	//decide on alt
	if(majorityIndex == referenceBase){
		calledGenotype = "0";

		//find second most common as alternative allele
		int second = pickIndexWithSecondHighestMetric(alleleCounts, 4, majorityIndex);
		altAlleles.push_back(second);
	} else {
		altAlleles.push_back(majorityIndex);
		calledGenotype = "1";
	}
};

void TCallerMajorityBase::callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//get per allele counts
	countAlleles(site);

	//now pick major among known alleles
	if(alleleCounts[referenceBase] > alleleCounts[altAlleles[0]]){
		calledGenotype = "0";
	} else if(alleleCounts[referenceBase] < alleleCounts[altAlleles[0]]){
		calledGenotype = "1";
	} else {
		//equal counts: pick at random
		if(randomGenerator->getRand() < 0.5)
			calledGenotype = "0";
		else
			calledGenotype = "1";
	}
};

/////////////////////////////////////////////////////////
// TCallerAllelePresence
/////////////////////////////////////////////////////////
TCallerAllelePresence::TCallerAllelePresence(TRandomGenerator* RandomGenerator):TCaller(RandomGenerator){
	//caller settings
	callerName = "Allele Presence Caller";
	filenameExtention = "_allelePresence.vcf";
	setAcceptableFields(&VCFGenotypeFields, "GT,DP,AD,GQ,AP");
	_usesPrior = true;

	//initialize allele counts
	MAP = -1;
};

void TCallerAllelePresence::fillPosteriors(TGenotypeLikelihoods & genotypeLikelihoods){
	//calculate posterior probabilities
	posterior.fill(genotypeLikelihoods, *genotypePrior);

	//sum for each base
	allelePostProb[0] = posterior[AA] + posterior[AC] + posterior[AG] + posterior[AT];
	allelePostProb[1] = posterior[AC] + posterior[CC] + posterior[CG] + posterior[CT];
	allelePostProb[2] = posterior[AG] + posterior[CG] + posterior[GG] + posterior[GT];
	allelePostProb[3] = posterior[AT] + posterior[CT] + posterior[GT] + posterior[TT];
};

void TCallerAllelePresence::callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(!priorSet) throw "Can not call AllelePresence genotypes: prior has not been set!";

	//fill posteriors for each allele
	fillPosteriors(genotypeLikelihoods);

	//find MAP
	MAP = pickIndexWithHighestMetric(allelePostProb, 4);

	//decide on alt
	if(MAP == referenceBase){
		calledGenotype = "0";

		//find second most common as alternative allele
		int second = pickIndexWithSecondHighestMetric(allelePostProb, 4, MAP);
		altAlleles.push_back(second);
	} else {
		altAlleles.push_back(MAP);
		calledGenotype = "1";
	}
};

void TCallerAllelePresence::callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(!priorSet) throw "Can not call AllelePresence genotypes: prior has not been set!";

	//fill posteriors for each allele
	fillPosteriors(genotypeLikelihoods);

	//find MAP
	double allelePostProbKnownAlleles[2];
	allelePostProbKnownAlleles[0] = allelePostProb[referenceBase];
	allelePostProbKnownAlleles[1] = allelePostProb[altAlleles[0]];

	int highest = pickIndexWithHighestMetric(allelePostProbKnownAlleles, 2);

	//decide on genotype (index 0 is ref base)
	if(highest == 0){
		calledGenotype = "0";
		MAP = referenceBase;
	} else {
		calledGenotype = "1";
		MAP = altAlleles[0];
	}
};

std::string TCallerAllelePresence::getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	return toString(qualMap.errorToPhredInt(1.0 - allelePostProb[MAP]));
};

std::string TCallerAllelePresence::getVCFGenotypeString_AP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	std::string ret = toString(qualMap.errorToPhredInt(posterior[0]));
	ret += ',' + toString(qualMap.errorToPhredInt(posterior[1]));
	ret += ',' + toString(qualMap.errorToPhredInt(posterior[2]));
	ret += ',' + toString(qualMap.errorToPhredInt(posterior[3]));
	return ret;
};

/////////////////////////////////////////////////////////
// TCallerDiploid
// common function between MLE, Bayes and gvcf callers
/////////////////////////////////////////////////////////
TCallerDiploid::TCallerDiploid(TRandomGenerator* RandomGenerator):TCaller(RandomGenerator){
	indexOfMax = -1;
	indexOfSecond = -1;
	imbalanceCalculated = false;
	missingGenotype = "./.";

	setAcceptableFields(&VCFGenotypeFields, "AB,AI");
};

void TCallerDiploid::clearAfterCall(){
	TCaller::clearAfterCall();
	imbalanceCalculated = false;
};

void TCallerDiploid::callGenotypeFromMetric(TGenotypeData & metric){
	indexOfMax = pickIndexWithHighestMetric(metric.pointerToData(), 10);
	indexOfSecond = pickIndexWithSecondHighestMetric(metric.pointerToData(), 10, indexOfMax);

	//decide on alternative alleles
	if(genoMap.genotypeToBase[indexOfMax][0] == referenceBase){
		if(genoMap.genotypeToBase[indexOfMax][1] == referenceBase){
			calledGenotype = "0/0";
			//MLE is homozygous reference -> find second best allele
			if(genoMap.genotypeToBase[indexOfSecond][0] == referenceBase)
				altAlleles.push_back(genoMap.genotypeToBase[indexOfSecond][1]);
			else if(genoMap.genotypeToBase[indexOfSecond][1] == referenceBase)
				altAlleles.push_back(genoMap.genotypeToBase[indexOfSecond][0]);
			else {
				//none is ref, pick at random
				int rand = randomGenerator->getRand() < 0.5 ? 0 : 1;
				altAlleles.push_back(genoMap.genotypeToBase[indexOfSecond][rand]);
			}
		} else {
			altAlleles.push_back(genoMap.genotypeToBase[indexOfMax][1]);
			calledGenotype = "0/1";
		}
	} else {
		if(genoMap.genotypeToBase[indexOfMax][1] == referenceBase){
			altAlleles.push_back(genoMap.genotypeToBase[indexOfMax][0]);
			calledGenotype = "0/1";
		} else {
			if(genoMap.genotypeToBase[indexOfMax][0] == genoMap.genotypeToBase[indexOfMax][1]){
				altAlleles.push_back(genoMap.genotypeToBase[indexOfMax][0]);
				calledGenotype = "1/1";

				//find second best allele, but give preference to reference in case likelihoods are equal
				if(_allowTriallelicSites){
					//int hetRef = genoMap.getGenotype(referenceBase, genoMap.genotypeToBase[indexOfMax][0]);
					//int homRef = genoMap.getGenotype(referenceBase, referenceBase);

					//only use second alternative allele in case het genotype with reference is less likely
					if(referenceBase == N || (metric[genoMap.getGenotype(referenceBase, genoMap.genotypeToBase[indexOfMax][0])] < metric[indexOfSecond] && metric[genoMap.getGenotype(referenceBase, referenceBase)] < metric[indexOfSecond])){
						if(genoMap.genotypeToBase[indexOfSecond][0] == referenceBase || genoMap.genotypeToBase[indexOfSecond][0] == altAlleles[0])
							altAlleles.push_back(genoMap.genotypeToBase[indexOfSecond][1]);
						else if(genoMap.genotypeToBase[indexOfSecond][1] == referenceBase || genoMap.genotypeToBase[indexOfSecond][1] == altAlleles[0])
							altAlleles.push_back(genoMap.genotypeToBase[indexOfSecond][0]);
						else {
							//none is ref, pick at random
							int rand = randomGenerator->getRand() < 0.5 ? 0 : 1;
							altAlleles.push_back(genoMap.genotypeToBase[indexOfSecond][rand]);
						}
					}
				}
			} else {
				if(_allowTriallelicSites){
					//allow triallelic sites
					altAlleles.push_back(genoMap.genotypeToBase[indexOfMax][0]);
					altAlleles.push_back(genoMap.genotypeToBase[indexOfMax][1]);
					calledGenotype = "1/2";
				} else {
					//decide on which of the two alternative alleles to pick -> check second highest
					if(genoMap.genotypeToBase[indexOfSecond][0] == genoMap.genotypeToBase[indexOfSecond][1]){
						if(genoMap.genotypeToBase[indexOfSecond][0] == genoMap.genotypeToBase[indexOfMax][0])
							altAlleles.push_back(genoMap.genotypeToBase[indexOfMax][0]);
						else if(genoMap.genotypeToBase[indexOfSecond][0] == genoMap.genotypeToBase[indexOfMax][1])
							altAlleles.push_back(genoMap.genotypeToBase[indexOfMax][0]);
						else {
							//neither alt 0 nor 1, pick at random
							int rand = randomGenerator->getRand() < 0.5 ? 0 : 1;
							altAlleles.push_back(genoMap.genotypeToBase[indexOfMax][rand]);
						}
					} else {
						if(genoMap.genotypeToBase[indexOfSecond][0] == referenceBase && (genoMap.genotypeToBase[indexOfSecond][1] == genoMap.genotypeToBase[indexOfMax][0] || genoMap.genotypeToBase[indexOfSecond][1] == genoMap.genotypeToBase[indexOfMax][1]))
							altAlleles.push_back(genoMap.genotypeToBase[indexOfSecond][1]);
						else if(genoMap.genotypeToBase[indexOfSecond][1] == referenceBase && (genoMap.genotypeToBase[indexOfSecond][0] == genoMap.genotypeToBase[indexOfMax][0] || genoMap.genotypeToBase[indexOfSecond][0] == genoMap.genotypeToBase[indexOfMax][1]))
							altAlleles.push_back(genoMap.genotypeToBase[indexOfSecond][0]);
						else {
							//neither alt 0 nor 1, pick at random
							int rand = randomGenerator->getRand() < 0.5 ? 0 : 1;
							altAlleles.push_back(genoMap.genotypeToBase[indexOfMax][rand]);
						}
					}

					//now call genotype from these alleles
					std::vector<int> indeces;
					callGenotypeFromMetricKnownAlleles(metric, indeces);
				}
			}
		}
	}
};


void TCallerDiploid::callGenotypeFromMetricKnownAlleles(const TGenotypeData & metric, std::vector<int> & indeces){
	//get genotypes
	int homRef = genoMap.getGenotype(referenceBase, referenceBase);
	int het = genoMap.getGenotype(referenceBase, altAlleles[0]);
	int homAlt = genoMap.getGenotype(altAlleles[0], altAlleles[0]);

	//find max
	double max = metric.at(homRef);
	if(metric.at(het) > max) max = metric.at(het);
	if(metric.at(homAlt) > max) max = metric.at(homAlt);

	//fill vector of all
	std::vector<std::string> vec;
	if(metric.at(homRef) == max){
		vec.push_back("0/0");
	}
	if(metric.at(het) == max){
		vec.push_back("0/1");
	}
	if(metric.at(homAlt) == max){
		vec.push_back("1/1");
	}

	calledGenotype = vec[randomGenerator->pickOne(vec.size())];
};

void TCallerDiploid::callGenotypeFromMetricKnownAllelesUpdateIndex(const TGenotypeData & metric){
	//initialize
	std::vector<std::string> gt;
	gt.push_back("0/0");
	gt.push_back("0/1");
	gt.push_back("1/1");

	//get genotypes
	int homRef = genoMap.getGenotype(referenceBase, referenceBase);
	int het = genoMap.getGenotype(referenceBase, altAlleles[0]);
	int homAlt = genoMap.getGenotype(altAlleles[0], altAlleles[0]);

	int indecesKnownAlleleGenotypes[3];
	indecesKnownAlleleGenotypes[0] = homRef;
	indecesKnownAlleleGenotypes[1] = het;
	indecesKnownAlleleGenotypes[2] = homAlt;

	double metricKnownAlleles[3];
	metricKnownAlleles[0] = metric.at(homRef);
	metricKnownAlleles[1] = metric.at(het);
	metricKnownAlleles[2] = metric.at(homAlt);

	int best = pickIndexWithHighestMetric(metricKnownAlleles, 3);
	int secondBest = pickIndexWithSecondHighestMetric(metricKnownAlleles, 3, indexOfMax);

	indexOfMax = indecesKnownAlleleGenotypes[best];
	calledGenotype = gt[best];
	indexOfSecond = indecesKnownAlleleGenotypes[secondBest];
};

std::string TCallerDiploid::getPerGenotypeMetricString(TGenotypeData & metric){
	//if you have alleles R, A, B, C then the order of the PL is: RR, RA, AA | RB, AB, BB | RC, AC, BC, CC
	//plot missing value (.) for all metrics involving the reference if the reference is N

	std::string ret;
	//first for reference base
	if(referenceBase == N)
		ret = ".";
	else
		ret = toString((int) metric[genoMap.genotypeMap[referenceBase][referenceBase]]);

	//now for alternative alleles
	if(altAlleles.size() > 0){
		if(referenceBase == N)
			ret += ",.";
		else
			ret += ',' + toString(metric[genoMap.genotypeMap[referenceBase][altAlleles[0]]]);
		ret += ',' + toString((int) metric[genoMap.genotypeMap[altAlleles[0]][altAlleles[0]]]);

		if(altAlleles.size() > 1){
			if(referenceBase == N)
				ret += ",.";
			else
				ret += ',' + toString(metric[genoMap.genotypeMap[referenceBase][altAlleles[1]]]);
			ret += ',' + toString(metric[genoMap.genotypeMap[altAlleles[0]][altAlleles[1]]]);
			ret += ',' + toString(metric[genoMap.genotypeMap[altAlleles[1]][altAlleles[1]]]);
		}

		if(altAlleles.size() > 2){
			if(referenceBase == N)
				ret += ",.";
			else
				ret += ',' + toString(metric[genoMap.genotypeMap[referenceBase][altAlleles[2]]]);
			ret += ',' + toString(metric[genoMap.genotypeMap[altAlleles[0]][altAlleles[2]]]);
			ret += ',' + toString(metric[genoMap.genotypeMap[altAlleles[1]][altAlleles[2]]]);
			ret += ',' + toString(metric[genoMap.genotypeMap[altAlleles[2]][altAlleles[2]]]);
		}
	}
	return ret;
};

void TCallerDiploid::calculateImbalance(const TSite & site){
	if(!imbalanceCalculated){
		if(!altAlleles.empty()){
			countAlleles(site);

			if(altAlleles.size() == 1){
				double sum = (alleleCounts[referenceBase] + alleleCounts[altAlleles[0]]);
				if(referenceBase == N || sum == 0){
					AB = '.'; AI = '.';
				} else {
					AB = toString(alleleCounts[referenceBase] / sum);
					AI = toString(randomGenerator->binomPValue(alleleCounts[referenceBase], alleleCounts[altAlleles[0]]));
				}
			} else {
				if(genoMap.genotypeToBase[indexOfMax][0] != genoMap.genotypeToBase[indexOfMax][1]){ //is het
					double sum = (double) alleleCounts[genoMap.genotypeToBase[indexOfMax][0]] + alleleCounts[genoMap.genotypeToBase[indexOfMax][1]];
					if(sum == 0){
						AB = '.'; AI = '.';
					} else {
						AB = toString(alleleCounts[genoMap.genotypeToBase[indexOfMax][0]] / sum);
						AI = toString(randomGenerator->binomPValue(alleleCounts[genoMap.genotypeToBase[indexOfMax][0]], alleleCounts[genoMap.genotypeToBase[indexOfMax][1]]));
					}
				} else { // is homo -> do it against the second alternative allele
					double sum = (double) alleleCounts[altAlleles[0]] + alleleCounts[altAlleles[1]];
					if(sum == 0){
						AB = '.'; AI = '.';
					} else {
						AB = toString(alleleCounts[altAlleles[0]] / sum);
						AI = toString(randomGenerator->binomPValue(alleleCounts[altAlleles[0]], alleleCounts[altAlleles[1]]));
					}
				}
			}
		}
		imbalanceCalculated = true;
	}
};

std::string TCallerDiploid::getVCFGenotypeString_AB(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(altAlleles.empty()) return ".";

	calculateImbalance(site);
	return AB;
};

std::string TCallerDiploid::getVCFGenotypeString_AI(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(altAlleles.empty()) return ".";

	calculateImbalance(site);
	return AI;
};

/////////////////////////////////////////////////////////
// TCallerMLE
/////////////////////////////////////////////////////////
TCallerMLE::TCallerMLE(TRandomGenerator* RandomGenerator):TCallerDiploid(RandomGenerator){
	//caller settings
	callerName = "MLE Caller";
	filenameExtention = "_MaximumLikelihood.vcf";
	setAcceptableFields(&VCFGenotypeFields, "AD,GQ,GL,PL");

	//set default tags to print
	printGenotypeFields("GT,DP,AD,GQ,PL");
};

void TCallerMLE::callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	callGenotypeFromMetric(genotypeLikelihoods);
};

void TCallerMLE::callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	callGenotypeFromMetricKnownAllelesUpdateIndex(genotypeLikelihoods);
};

std::string TCallerMLE::getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	return toString(qualMap.errorToPhredInt(genotypeLikelihoods[indexOfSecond] / genotypeLikelihoods[indexOfMax]));
};

std::string TCallerMLE::getVCFGenotypeString_GL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//normalize
	for(uint8_t g=0; g<10; ++g)
		tmpGenoData[g] = log10(genotypeLikelihoods[g] / genotypeLikelihoods[indexOfMax]);

	//get string
	return getPerGenotypeMetricString(tmpGenoData);
};

std::string TCallerMLE::getVCFGenotypeString_PL(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//normalize
	double phredMax = qualMap.errorToPhred(genotypeLikelihoods[indexOfMax]);
	for(uint8_t g=0; g<10; ++g)
		tmpGenoData[g] = (int) round(qualMap.errorToPhred(genotypeLikelihoods[g]) - phredMax);

	//get string
	return getPerGenotypeMetricString(tmpGenoData);
};

//------------------------------------------------------
// TCallerBayes
//------------------------------------------------------
TCallerBayes::TCallerBayes(TRandomGenerator* RandomGenerator):TCallerDiploid(RandomGenerator){
	//caller settings
	callerName = "Bayesian Caller";
	filenameExtention = "_MaximumAPosteriori.vcf";
	setAcceptableFields(&VCFGenotypeFields, "AD,GQ,GP");
	printGenotypeFields("GT,DP,AD,GQ,GP");
	_usesPrior = true;
};


void TCallerBayes::callGenotype(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(!priorSet) throw "Can not call Bayesian genotypes: prior has not been set!";

	//calculate posterior probabilities
	posterior.fill(genotypeLikelihoods, *genotypePrior);

	//call
	callGenotypeFromMetric(posterior);
};

void TCallerBayes::callGenotypeKnownAlleles(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	if(!priorSet) throw "Can not call Bayesian genotypes: prior has not been set!";

	//calculate posterior probabilities
	posterior.fill(genotypeLikelihoods, *genotypePrior);

	//call
	callGenotypeFromMetricKnownAllelesUpdateIndex(posterior);
};

std::string TCallerBayes::getVCFGenotypeString_GQ(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	return toString(qualMap.errorToPhredInt(1.0 - posterior[indexOfMax]));
};

std::string TCallerBayes::getVCFGenotypeString_GP(const TSite & site, TGenotypeLikelihoods & genotypeLikelihoods){
	//phred
	for(uint8_t g=0; g<10; ++g)
		tmpGenoData[g] = (int) qualMap.errorToPhredInt(posterior[g]);

	//get string
	return getPerGenotypeMetricString(tmpGenoData);
};

