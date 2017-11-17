/*
 * TAlignmentParser.cpp
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#include "TAlignmentParser.h"

TAlignmentParser::TAlignmentParser(){
	maxSize = 0;
	readGroupTable = NULL;
	logfile = NULL;
	_keepDuplicates = false;
	initialized = false;

	//details
	length = 0;
	chrNumber = -1;
	readGroupId = -1;
	position = 0;
	isReverseStrand = false;
	passedFilters = false;
	parsed = false;
	changed = false;

	//per base data
	base = NULL;
	baseAsChar = NULL;
	context = NULL;
	quality = NULL;
	qualityOriginal = NULL;
	errorRates = NULL;
	recalibrated = false;
	qualityRecalibrated = NULL;
	aligned = NULL;
	alignedPos = NULL;
	distFrom3Prime = NULL;
	distFrom5Prime = NULL;

	//soft clipped data
	softClippedEntry = 0;
	softClippedLength = NULL;
	softClippedBase = NULL;
	softClippedQuality = NULL;

	//tmp variables
	i = 0;
	d = 0;
	k = 0;
	p = 0;
}

TAlignmentParser::TAlignmentParser(TReadGroups* ReadGroupTable, unsigned int MaxSize, TLog* Logfile){
	TAlignmentParser();
	init(ReadGroupTable, MaxSize, Logfile);
};

void TAlignmentParser::init(TReadGroups* ReadGroupTable, unsigned int MaxSize, TLog* Logfile){
	clear();

	maxSize = MaxSize;
	readGroupTable = ReadGroupTable;
	logfile = Logfile;

	//data
	base = new Base[maxSize];
	baseAsChar = new char[maxSize];
	context = new BaseContext[maxSize];
	qualityOriginal = new int[maxSize];
	qualityRecalibrated = new int[maxSize];
	errorRates = new double[maxSize];
	aligned = new bool[maxSize];
	alignedPos = new int[maxSize];
	distFrom3Prime = new int[maxSize];
	distFrom5Prime = new int[maxSize];

	//soft clipped data
	softClippedLength = new int[2];
	softClippedBase = new char*[2];
	softClippedQuality = new char*[2];
	softClippedBase[0] = new char[maxSize];
	softClippedBase[1] = new char[maxSize];
	softClippedQuality[0] = new char[maxSize];
	softClippedQuality[1] = new char[maxSize];

	initialized = true;
};

void TAlignmentParser::clear(){
	if(initialized){
		delete[] base;
		delete[] qualityOriginal;
		delete[] qualityRecalibrated;
		delete[] errorRates;
		delete[] aligned;
		delete[] alignedPos;
		delete[] distFrom3Prime;
		delete[] distFrom5Prime;

		delete[] softClippedLength;
		delete[] softClippedBase[0];
		delete[] softClippedBase[1];
		delete[] softClippedBase;
		delete[] softClippedQuality[0];
		delete[] softClippedQuality[1];
		delete[] softClippedQuality;

		initialized = false;
	}
};

void TAlignmentParser::parseBasesQualities(){
	// iterate over CigarOps
	d = 0; //index regarding data structures
	k = 0; //index inside read
	p = 0; //index regarding reference position (!= k for indels)
	softClippedEntry = 0; //softclipped bases to be added
	softClippedLength[0] = 0; softClippedLength[1] = 0;

	std::vector<BamTools::CigarOp>::const_iterator cigarIter = bamAlignment.CigarData.begin();
	std::vector<BamTools::CigarOp>::const_iterator cigarEnd  = bamAlignment.CigarData.end();

	for ( ; cigarIter != cigarEnd; ++cigarIter ) {
		const BamTools::CigarOp& op = (*cigarIter);
		switch ( op.Type ) {

			// for 'M', '=' or 'X': just copy
			case (BamTools::Constants::BAM_CIGAR_MATCH_CHAR)    :
			case (BamTools::Constants::BAM_CIGAR_SEQMATCH_CHAR) :
			case (BamTools::Constants::BAM_CIGAR_MISMATCH_CHAR) :
				for(i=0; i<op.Length; ++i, ++d, ++k, ++p){
					base[d] = genoMap.getBaseOnlyCapitals(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					qualityOriginal[d] = (int) bamAlignment.Qualities[k];
					errorRates[d] = qualityMap.charToErrorMap[bamAlignment.Qualities[k]];
					aligned[d] = true;
					alignedPos[d] = p;
				}
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
				break;

			//for 'S' - soft clip: ignore bases, but increase k
			case (BamTools::Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
				//add bases to softclipped entries
				for(i=0; i<op.Length; ++i, ++d, ++k, ++p){
					softClippedBase[softClippedEntry][softClippedLength[softClippedEntry]] = bamAlignment.QueryBases[k];
					softClippedQuality[softClippedEntry][softClippedLength[softClippedEntry]] = bamAlignment.Qualities[k];
					++softClippedLength[softClippedEntry];
				}
				break;

			//for 'I' - insertion: copy bases, but put aligned pos to
			case (BamTools::Constants::BAM_CIGAR_INS_CHAR)      :
				for(i=0; i<op.Length; ++i, ++d, ++k){
					base[d] = genoMap.getBaseOnlyCapitals(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					qualityOriginal[d] = (int) (char) bamAlignment.Qualities[k];
					errorRates[d] = qualityMap.charToErrorMap[bamAlignment.Qualities[k]];
					aligned[d] = false;
					alignedPos[d] = -1;
				}
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
				break;


			// for 'D' - deletion: just add to postion
			case (BamTools::Constants::BAM_CIGAR_DEL_CHAR) :
				p += op.Length;
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
				break;

			// for 'N' - skipped region: copy but say that bases were not aligned
			case (BamTools::Constants::BAM_CIGAR_REFSKIP_CHAR) :
				for(i=0; i<op.Length; ++i, ++d, ++k, ++p){
					base[d] = genoMap.getBaseOnlyCapitals(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					qualityOriginal[d] = (int) bamAlignment.Qualities[k];
					errorRates[d] = qualityMap.charToErrorMap[bamAlignment.Qualities[k]];
					aligned[d] = false;
					alignedPos[d] = p;
				}
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
				break;

			// for 'H' - hard clip: do nothing as these bases are not present in SEQ
			case (BamTools::Constants::BAM_CIGAR_HARDCLIP_CHAR) :
				break;

			// invalid CIGAR op-code
			default:
				throw (std::string) "CIGAR operation type '" + op.Type + "' not supported!";
		}
	}

	//update length
	length = k;

	if(length != bamAlignment.Length)
		throw "Length mismatch!";
};

void TAlignmentParser::setDistancesFromEnds(){
	//is it paired-end?
	if(bamAlignment.IsProperPair()){
		if(bamAlignment.IsReverseStrand()){
			//reverse (can be either first or second mate, but it's the one that comes second in bam file)
			//hence distance from 3' is given by f(dist since beginning of fragment) = f(insert - len + pos)
			//and distance from 5' is given as f(end of fragment) = f(len - pos - 1)
			k = abs(bamAlignment.InsertSize) - length;
			p = length - 1;
			for(d=0; d<length; ++d){
				distFrom3Prime[d] = k + d;
				distFrom5Prime[d] = p - d;
			}
		} else {
			//forward (can be either first or second mate, but it's the one that comes first in bam file)
			//Hence distance from 3' is given as a function of pos
			//And distance from 5' is given by (length of fragment) - pos -1
			//NOTE! we ignore indels when calculating distance from 5' since we can not know this info.
			//Luckily, this has only minimal effect since these distances are far from fragment ends
			p = abs(bamAlignment.InsertSize) - 1;
			for(d=0; d<length; ++d){
				distFrom3Prime[d] = d;
				distFrom5Prime[d] = p - d;
			}
		}
	} else {
		//treat as single end
		p = length - 1;
		if(bamAlignment.IsReverseStrand()){
			//not in pair & reverse
			//Hence distance from 5' is just pos
			//And distance from 3' is just len - pos - 1
			for(d=0; d<length; ++d){
				distFrom3Prime[d] = p - d;
				distFrom5Prime[d] = d;
			}

		} else {
			//not in pair & forward
			//Hence distance from 3' is given as a function of pos
			//And distance from 5' is given by len - pos - 1
			for(d=0; d<length; ++d){
				distFrom3Prime[d] = d;
				distFrom5Prime[d] = p - d;
			}
		}
	}
};

void TAlignmentParser::fillContext(){
	if(isReverseStrand){
		//reverse
		for(d=0; d<(length-1); ++d)
			context[d] = genoMap.contextMap[base[d+1]][base[d]];
		context[d] = genoMap.contextMap[N][base[d]];
	} else {
		//forward
		context[0] = genoMap.contextMap[N][base[d]];
		for(d=1; d<length; ++d)
			context[d] = genoMap.contextMap[base[d-1]][base[d]];
	}
};

//------------------------------
//public functions
//------------------------------
bool TAlignmentParser::readAlignment(BamTools::BamReader & bamReader){
	if(!bamReader.GetNextAlignment(bamAlignment))
		return false;

	//add basic info
	parsed = false;
	recalibrated = false;
	changed = false;
	quality = qualityOriginal;
	chrNumber = bamAlignment.RefID;
	position = bamAlignment.Position;
	isReverseStrand = bamAlignment.IsReverseStrand();

	//Extract Read Group Info
	bamAlignment.GetTag("RG", readGroup);
	readGroupId = readGroupTable->find(readGroup);

	//check if read passes basic QC
	passedFilters = bamAlignment.IsMapped() && !bamAlignment.IsFailedQC() && bamAlignment.IsPrimaryAlignment() && (_keepDuplicates || !bamAlignment.IsDuplicate());

	//check read length
	if(bamAlignment.AlignedBases.size() > maxSize)
		throw "Alignment '" +  bamAlignment.Name + "' is longer than the max read length! Please change max read length to parse this data.";

	//check if insert size is shorter than read, this means we are reading the adaptor sequence
	if(bamAlignment.IsPaired() && abs(bamAlignment.InsertSize) <= bamAlignment.AlignedBases.length()){
		logfile->warning("The following alignment is longer than its insert size: " + bamAlignment.Name);
		passedFilters = false;
	}

	return true;
}


void TAlignmentParser::parse(){
	if(!parsed){
		//first parse bases and qualities
		parseBasesQualities();

		//then update distances from ends
		setDistancesFromEnds();

		parsed = true;
	}
};

void TAlignmentParser::recalibrate(TRecalibration & recalObject){
	//make sure read is parsed
	parse();

	if(recalObject.recalibrationChangesQualities()){
		//recalibrate quality scores
		for(d=0; d<length; ++d){
			k = length - d - 1;
			errorRates[d] = recalObject.getErrorRate(readGroupId, qualityOriginal[d], d, k, context[d]);
			qualityRecalibrated[d] = qualityMap.errorToPhred(errorRates[d]);
		}
		quality = qualityRecalibrated;
		changed = true;
	} else changed = false;
	recalibrated = true;

};

void TAlignmentParser::recalibrate(TRecalibration & recalObject, TPMD* pmdObjects, BamTools::Fasta & reference){
	//make sure read is parsed
	parse();

	//get reference sequence
	reference.GetSequence(chrNumber, position, alignedPos[length-1], tmpString);

	//recalibrate quality scores
	for(d=0; d<length; ++d){
		k = length - d - 1;
		if(recalObject.recalibrationChangesQualities())
			errorRates[d] = recalObject.getErrorRate(readGroupId, qualityOriginal[d], d, k, context[d]);

		//now add effect of PMD
		if(aligned[d]){
			if(base[d] == T && tmpString[alignedPos[d]] == 'C')
				errorRates[d] = 1.0 - ((1.0 - errorRates[d])*(1.0 - pmdObjects[readGroupId].getProbCT(distFrom5Prime[d]))); //this is mapDamage2, Krishna: qual*(1-pmdCT) + (1-qual)*pmdCT;
			else if(base[d] == A && tmpString[alignedPos[d]] == 'G')
				errorRates[d] = 1.0 - ((1.0 - errorRates[d])*(1.0 - pmdObjects[readGroupId].getProbGA(distFrom3Prime[d]))); //this is mapDamage2, Krishna: qual*(1-pmdGA) + (1-qual)*pmdGA;
		}

		qualityRecalibrated[d] = qualityMap.errorToPhred(errorRates[d]);
	}

	//set pointer to recalibrated scores
	quality = qualityRecalibrated;
	recalibrated = true;
	changed = true;
};

void TAlignmentParser::save(BamTools::BamWriter & bamWriter){
	if(!changed){
		if(!bamWriter.SaveAlignment(bamAlignment))
			throw "Read '" + bamAlignment.Name + "' could not be written!";
	}
	else {
		//means that read has been modified.
		//Currently quality recalibration is the only possible change.
		//But will need to think how to deal with merging and such...

		//assume that only bases and quality scores where changed
		tmpString.clear();
		tmpString2.clear();
		for(d=0; d<length; ++d){
			tmpString += genoMap.baseToChar[base[d]];
			tmpString2 += (char) quality[d];
		}

		bamAlignment.QueryBases = tmpString;
		bamAlignment.Qualities = tmpString2;

		if(!bamWriter.SaveAlignment(bamAlignment))
			throw "Read '" + bamAlignment.Name + "' could not be written!";
	}
};

void TAlignmentParser::print(){
	std::cout << "NAME:\t" << bamAlignment.Name << std::endl;
	std::cout << "LEN:\t" << length << std::endl;

	//print bases
	std::cout << "SEQ:\t";
	for(d=0; d<length; ++d)
		std::cout << genoMap.getBaseAsChar(base[d]);
	std::cout << std::endl;

	//print qualities
	std::cout << "QUAL:\t";
	for(d=0; d<length; ++d)
		std::cout << (char) (quality[d]+33);
	std::cout << std::endl;

	//print aligned pos
	std::cout << "POS:\t";
	for(d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		if(aligned[d])
			std::cout << alignedPos[d];
		else
			std::cout << "-";
	}
	std::cout << std::endl;

	//print dist from 3'
	std::cout << "dist 3':\t";
	for(d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		std::cout << distFrom3Prime[d];
	}
	std::cout << std::endl;

	//print dist from 5'
	std::cout << "dist 5':\t";
	for(d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		std::cout << distFrom5Prime[d];
	}
	std::cout << std::endl;
}




