/*
 * TAlignmentParser.cpp
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#include "TAlignmentParser.h"

//-----------------------------------------------------
//TFastaBuffer
//-----------------------------------------------------
TFastaBuffer::TFastaBuffer(BamTools::Fasta* Reference){
	bufferSize = 10000000;
	reference = Reference;
	referenceSequence = "";
	curStart = -1;
	curChr = -1;
	curEnd = -1;
}


void TFastaBuffer::moveTo(const int & chr, const int32_t & pos){
	curChr = chr;
	curStart = pos;
	curEnd = pos + bufferSize;
	if(!reference->GetSequence(chr, curStart, curEnd, referenceSequence))
		throw "Problem reading " + toString(chr) + ":" + toString(curStart) + "-" + toString(curEnd) + " from fasta file!";
}

void TFastaBuffer::fill(const int & chr, const int32_t & start, const int32_t end, std::string & ref){
	//move buffer, if necessary
	if(chr != curChr || end > curEnd || start < curStart)
		moveTo(chr, start);

	//now copy to string
	ref.assign(referenceSequence, start - curStart, end - start + 1);
}

//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
TAlignmentParser::TAlignmentParser(){
	maxSize = 0;
	readGroupTable = NULL;
	logfile = NULL;
	_keepDuplicates = false;
	initialized = false;
	applyQualityFilter = false;
	minQual = 33;
	maxQual = 126;
	minQualForPrinting = 33;
	maxQualForPrinting = 126;

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
	pmdCT = NULL;
	pmdGA = NULL;

	//soft clipped data
	softClippedEntry = 0;
	softClippedLength = NULL;
	softClippedBase = NULL;
	softClippedQuality = NULL;

	//reference
	hasReference = false;
	fastaBuffer = NULL;

	//tmp variables
	i = 0;
	d = 0;
	k = 0;
	p = 0;
};

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
	pmdCT = new double[maxSize];
	pmdGA = new double[maxSize];

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
		delete[] baseAsChar;
		delete[] context;
		delete[] qualityOriginal;
		delete[] qualityRecalibrated;
		delete[] errorRates;
		delete[] aligned;
		delete[] alignedPos;
		delete[] distFrom3Prime;
		delete[] distFrom5Prime;
		delete[] pmdCT;
		delete[] pmdGA;

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

void TAlignmentParser::setQualityFilters(int MinQual, int MaxQual){
	applyQualityFilter = true;
	minQual = MinQual;
	maxQual = MaxQual;
};

void TAlignmentParser::setQualityRangeForPrinting(int minQual, int maxQual){
	minQualForPrinting = minQual;
	maxQualForPrinting = maxQual;
};

void TAlignmentParser::parseBasesQualities(){
	// iterate over CigarOps
	d = 0; //index regarding data structures
	k = 0; //index inside read
	p = 0; //index regarding reference position (!= k for indels)
	softClippedEntry = 0; //softclipped bases to be added
	softClippedLength[0] = 0; softClippedLength[1] = 0;

	cigarIter = bamAlignment.CigarData.begin();
	cigarEnd  = bamAlignment.CigarData.end();
	int counter = 0;
	for( ; cigarIter != cigarEnd; ++cigarIter, ++counter ){
		const BamTools::CigarOp& op = (*cigarIter);
		switch ( op.Type ) {

			// for 'M', '=' or 'X': just copy
			case (BamTools::Constants::BAM_CIGAR_MATCH_CHAR)    :
			case (BamTools::Constants::BAM_CIGAR_SEQMATCH_CHAR) :
			case (BamTools::Constants::BAM_CIGAR_MISMATCH_CHAR) :
				for(i=0; i<op.Length; ++i, ++d, ++k, ++p){
					base[d] = genoMap.getBase(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					qualityOriginal[d] = (int) bamAlignment.Qualities[k];
					errorRates[d] = qualityMap.qualityToErrorMap[(int) bamAlignment.Qualities[k]];
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
					//need to initialize quality for quality filter and bases for context
					base[d] = genoMap.getBase(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					qualityOriginal[d] = -1;
					errorRates[d] = qualityMap.qualityToErrorMap[(int) bamAlignment.Qualities[k]];
					aligned[d] = false;
					alignedPos[d] = -1;
				}
				break;

			//for 'I' - insertion: copy bases, but put aligned pos to
			case (BamTools::Constants::BAM_CIGAR_INS_CHAR)      :
				for(i=0; i<op.Length; ++i, ++d, ++k){
					base[d] = genoMap.getBase(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					qualityOriginal[d] = (int) (char) bamAlignment.Qualities[k];
					errorRates[d] = qualityMap.qualityToErrorMap[(int) bamAlignment.Qualities[k]];
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
					base[d] = genoMap.getBase(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					qualityOriginal[d] = (int) bamAlignment.Qualities[k];
					errorRates[d] = qualityMap.qualityToErrorMap[(int) bamAlignment.Qualities[k]];
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

void TAlignmentParser::filterForBaseQuality(){
	//set base to N if outside quality filter
	for(d=0; d<length; ++d){
		if(qualityOriginal[d] < minQual || qualityOriginal[d] > maxQual){
			base[d] = N;
		}
	}
};

void TAlignmentParser::filterForPrintingBaseQuality(std::string & qual){
	//set base to N if outside quality filter
	for(stringIt = qual.begin() ; stringIt < qual.end(); ++stringIt){

		if((int) *stringIt < minQualForPrinting)
			*stringIt = (char) minQualForPrinting;
		else if((int) *stringIt > maxQualForPrinting)
			*stringIt = (char) maxQualForPrinting;
	}
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
				distFrom5Prime[d] = p - d;
				distFrom3Prime[d] = k + d;
			}
		} else {
			//forward (can be either first or second mate, but it's the one that comes first in bam file)
			//Hence distance from 3' is given as a function of pos
			//And distance from 5' is given by (length of fragment) - pos -1
			//NOTE! we ignore indels when calculating distance from 5' since we can not know this info.
			//Luckily, this has only minimal effect since these distances are far from fragment ends
			p = abs(bamAlignment.InsertSize) - 1;
			for(d=0; d<length; ++d){
				distFrom5Prime[d] = d;
				distFrom3Prime[d] = p - d;
			}
		}
	} else {
		//treat as single end
		p = length - 1;
		if(bamAlignment.IsReverseStrand()){
			//not in pair & reverse
			//Hence distance from 3' is just pos
			//And distance from 5' is just len - pos - 1
			for(d=0; d<length; ++d){
				distFrom5Prime[d] = p - d;
				distFrom3Prime[d] = d;
			}

		} else {
			//not in pair & forward
			//Hence distance from 5' is given as a function of pos
			//And distance from 3' is given by len - pos - 1
			for(d=0; d<length; ++d){
				distFrom5Prime[d] = d;
				distFrom3Prime[d] = p - d;
			}
		}
	}
};

void TAlignmentParser::fillContext(){
	if(isReverseStrand){
		//reverse
		for(d=0; d<(length-1); ++d){
//			std::cout << "getting Context for " << base[d-1] << ", " << base[d] << std::flush;
			//std::cout << " -> " << genoMap.contextMap[base[d-1]][base[d]] << std::endl;
			context[d] = genoMap.contextMap[base[d+1]][base[d]];
		}
		context[d] = genoMap.contextMap[N][base[d]];
	} else {
		//forward
		context[0] = genoMap.contextMap[N][base[0]];
		for(d=1; d<length; ++d){
//			std::cout << "getting Context for " << base[d-1] << ", " << base[d] << std::flush;
			//std::cout << " -> " << genoMap.contextMap[base[d-1]][base[d]] << std::endl;
			context[d] = genoMap.contextMap[base[d-1]][base[d]];
		}
	}
};

void TAlignmentParser::fillPmdProbabilities(TPMD* pmdObjects){
	for(d=0; d<length; ++d){
		pmdCT[d] = pmdObjects[readGroupId].getProbCT(distFrom5Prime[d]);
		pmdGA[d] = pmdObjects[readGroupId].getProbGA(distFrom3Prime[d]);
	}
};

void TAlignmentParser::addReference(BamTools::Fasta* reference){
	hasReference = true;
	fastaBuffer = new TFastaBuffer(reference);
};

void TAlignmentParser::fillReferenceSequence(){
	if(!hasReference) //is this check really necessary?
		throw "No reference provided!";

	fastaBuffer->fill(chrNumber, position, position + alignedPos[length-1], referenceSequence);
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
	isProperPair = bamAlignment.IsProperPair();

	//Extract Read Group Info
	bamAlignment.GetTag("RG", readGroup);
	//TODO: add check of whether RG is used
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

		//fill context for each base
		fillContext();

		//apply filters
		//TODO: should that be on the recalibrated quality scores instead???
		if(applyQualityFilter)
			filterForBaseQuality();

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
			qualityRecalibrated[d] = qualityMap.errorToQuality(errorRates[d]);
		}

		quality = qualityRecalibrated;
		changed = true;
	} else changed = false;
	recalibrated = true;

};

void TAlignmentParser::recalibrate(TRecalibration & recalObject, TPMD* pmdObjects){
	//make sure read is parsed
	parse();

	//get reference sequence
	fillReferenceSequence();

	//get PMD probs
	fillPmdProbabilities(pmdObjects);

	//recalibrate quality scores
	for(d=0; d<length; ++d){
		k = length - d - 1;
		if(recalObject.recalibrationChangesQualities())
			errorRates[d] = recalObject.getErrorRate(readGroupId, qualityOriginal[d], d, k, context[d]);

		//now add effect of PMD
		if(aligned[d]){
			if(base[d] == T && referenceSequence[alignedPos[d]] == 'C')
				errorRates[d] = 1.0 - ((1.0 - errorRates[d])*(1.0 - pmdCT[d])); //this is mapDamage2, Krishna: qual*(1-pmdCT) + (1-qual)*pmdCT;
			else if(base[d] == A && referenceSequence[alignedPos[d]] == 'G')
				errorRates[d] = 1.0 - ((1.0 - errorRates[d])*(1.0 - pmdGA[d])); //this is mapDamage2, Krishna: qual*(1-pmdGA) + (1-qual)*pmdGA;
		}

		qualityRecalibrated[d] = qualityMap.errorToQuality(errorRates[d]);
	}

	//set pointer to recalibrated scores
	quality = qualityRecalibrated;
	recalibrated = true;
	changed = true;
};

void TAlignmentParser::binQualityScores(){
	//make sure read is parsed
	parse();

	//bin quality scores as done by Illumina
	for(d=0; d<length; ++d){
		quality[d] = qualityMap.illuminaQualityBins[quality[d]];
	}
	changed = true;
};


void TAlignmentParser::addToPMDTables(TPMDTables & pmdTables){
	//make sure read is parsed
	parse();

	//get reference sequence
	fillReferenceSequence();

	//tmp variables
	Base ref, read;

	//check if it is forward or reverse strand!
	if(isReverseStrand){
		for(d=0; d<length; ++d){
			if(aligned[d] && base[d] != N){
				ref = genoMap.flipBase(referenceSequence[alignedPos[d]]);
				read = genoMap.baseToFlippedBase[base[d]];
				pmdTables.addForward(readGroupId, distFrom3Prime[d], ref, read);
				pmdTables.addReverse(readGroupId, distFrom5Prime[d], ref, read);
			}
		}
	} else {
		for(d=0; d<length; ++d){
			if(aligned[d] && base[d] != N){
				ref = genoMap.getBase(referenceSequence[alignedPos[d]]);
				pmdTables.addForward(readGroupId, distFrom5Prime[d], ref, base[d]);
				pmdTables.addReverse(readGroupId, distFrom3Prime[d], ref, base[d]);
			}
		}
	}
};

double TAlignmentParser::calculatePMDS(double & pi, TPMD* pmdObjects){
	//make sure read is parsed
	parse();

	//variables
	double PMDS = 0.0;
	double probPMD, probNoPMD;
	double epsThird;
	double fourEpsThird;

	//get reference
	fillReferenceSequence();

	//get PMD probs
	fillPmdProbabilities(pmdObjects);

	//go over all bases in read
	for(d=0; d<length; ++d){
		//limit to aligned positions
		if(aligned[d]){
			//Prepare variables
			epsThird = errorRates[d] / 3.0;
			fourEpsThird = 4.0 * epsThird;

			//calc likelihoods
			if(referenceSequence[alignedPos[d]] == 'A'){
				if(base[d] == A){
					probPMD = 1.0 - errorRates[d] - pi + fourEpsThird*pi + pmdGA[d]*pi/3.0*(1.0-fourEpsThird);
					probNoPMD = 1.0 - errorRates[d] - pi + fourEpsThird*pi;
				}
				else if(base[d] == C){ //ok
					probPMD = errorRates[d] - fourEpsThird*pi + pi - pi*pmdCT[d]*(fourEpsThird-1.0);
					probNoPMD = errorRates[d] - fourEpsThird*pi + pi;
				}
				else if(base[d] == G){
					probPMD = errorRates[d] - fourEpsThird*pi + pi + pi*pmdGA[d]*(fourEpsThird-1.0);
					probNoPMD = errorRates[d] - fourEpsThird*pi + pi;
				}
				else if(base[d] == T){
					probPMD = errorRates[d] - fourEpsThird*pi + pi + pi*pmdCT[d]*(1.0-fourEpsThird);
					probNoPMD = errorRates[d] - fourEpsThird*pi + pi;
				}
			}
			else if (referenceSequence[alignedPos[d]] == 'C'){
				if(base[d] == A){
					probPMD = errorRates[d] + pi - 2.0*errorRates[d]*pi + pi*pmdGA[d]*(1.0-fourEpsThird);
					probNoPMD = errorRates[d] + pi - 2.0*errorRates[d]*pi;
				}
				else if(base[d] == C){
					probPMD = 1.0 - pi - errorRates[d] + fourEpsThird*pi + (1.0-pi)*pmdCT[d]*(fourEpsThird-1.0);
					probNoPMD = 1.0 - pi - errorRates[d] + fourEpsThird*pi;
				}
				else if(base[d] == G){
					probPMD = errorRates[d] - fourEpsThird*pi + pi + pmdGA[d]*(fourEpsThird*pi - pi);
					probNoPMD = errorRates[d] - fourEpsThird*pi + pi;
				}
				else if(base[d] == T){
					probPMD = epsThird + (1.0-pi)*pmdCT[d]*(1.0-fourEpsThird);
					probNoPMD = epsThird;
				}
			}
			else if (referenceSequence[alignedPos[d]] == 'G'){
				if(base[d] == A){
					probPMD = pmdGA[d]*(3.0-3.0*pi+4.0*errorRates[d]+4.0*errorRates[d]*pi) + errorRates[d] - fourEpsThird*pi + pi;
					probNoPMD = errorRates[d] - fourEpsThird*pi + pi;
				}
				else if(base[d] == C){
					probPMD = errorRates[d] - fourEpsThird*pi + pi + pi*pmdCT[d]*(fourEpsThird - 1.0);
					probNoPMD = errorRates[d] - fourEpsThird*pi + pi;
				}
				else if(base[d] == G){
					probPMD = 1.0 - errorRates[d] - pi + fourEpsThird*pi + (1.0-pi)*pmdGA[d]*(fourEpsThird-1.0);
					probNoPMD = 1.0 - errorRates[d] - pi + fourEpsThird*pi;
				}
				else if(base[d] == T){
					probPMD = errorRates[d] - fourEpsThird*pi + pi + pi*pmdCT[d]*(1.0-fourEpsThird);
					probNoPMD = errorRates[d] - fourEpsThird*pi + pi;
				}
			}
			else if(referenceSequence[alignedPos[d]] == 'T'){
				if(base[d] == A){
					probPMD = errorRates[d] - fourEpsThird*pi + pi - epsThird*pi*pmdCT[d] + pi*pmdGA[d]*(1.0-errorRates[d]);
					probNoPMD = errorRates[d] - fourEpsThird*pi + pi;
				}
				else if(base[d] == C){
					probPMD = errorRates[d] - fourEpsThird*pi + pi + pi*pmdCT[d]*(fourEpsThird - 1.0);
					probNoPMD = errorRates[d] - fourEpsThird*pi + pi;
				}
				else if(base[d] == G){
					probPMD = errorRates[d] - fourEpsThird*pi + pi + pi*pmdGA[d]*(fourEpsThird - 1.0);
					probNoPMD = errorRates[d] - fourEpsThird*pi + pi;
				}
				else if(base[d] == T){
					probPMD = 1.0 - errorRates[d] - pi + fourEpsThird*pi + pmdCT[d]*(pi/3.0-fourEpsThird*pi/3.0);
					probNoPMD = 1.0 - errorRates[d] - pi + fourEpsThird*pi;
				}
			}

			//now add to PMDS
			PMDS += log(probPMD/probNoPMD);
		}
	}

	return PMDS;
}

void TAlignmentParser::assessSoftClipping(int & S_left, int & middle, int & S_right){
	//count S, not S, S pattern from cigar string
	S_left = 0;
	S_right = 0;
	middle = 0;
	static bool reachedMiddle = false;

	cigarIter = bamAlignment.CigarData.begin();
	cigarEnd  = bamAlignment.CigarData.end();

	for( ; cigarIter != cigarEnd; ++cigarIter ){
		if(cigarIter->Type == 'S'){
			if(reachedMiddle) S_right += cigarIter->Length;
			else S_left += cigarIter->Length;
		} else {
			reachedMiddle = true;
			middle += cigarIter->Length;
		}
	}
};

void TAlignmentParser::addToQualityTable(TQualityTable & qualTable){
	parse();
	qualTable.add(quality, length);
};

//--------------------------------------------
//functions to modify alignment
//--------------------------------------------
void TAlignmentParser::updateOptionalSamField(std::string tag, float value){
	if(bamAlignment.HasTag(tag) == false) bamAlignment.AddTag(tag, "f", value);
	else bamAlignment.EditTag(tag, "f", value);
};

void TAlignmentParser::downsampleAlignment(double& fraction, TRandomGenerator& randomGenerator){
	for(int d=0; d<length; ++d){
		double r = randomGenerator.getRand();
		if(r < fraction){
			base[d] = N;
			quality[d] = 0;
		}
	}
	changed = true;
}
//--------------------------------------------
//functions to write / print alignment
//--------------------------------------------
void TAlignmentParser::save(BamTools::BamWriter & bamWriter){
	if(changed){
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
	}

	//make sure quality are printed within range
	filterForPrintingBaseQuality(bamAlignment.Qualities);

	//now write alignment
	if(!bamWriter.SaveAlignment(bamAlignment))
		throw "Read '" + bamAlignment.Name + "' could not be written!";
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
		std::cout << (char) quality[d];
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




