/*
 * TAlignment.cpp
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */

#include "TAlignment.h"

TAlignment::TAlignment(){
	//details
	empty = true;
	maxSize = 0;
	length = 0;
	chrNumber = -1;
	readGroupId = -1;
	position = 0;
	isReverseStrand = false;
	isProperPair = false;
	passedFilters = false;
	parsed = false;
	changed = false;
	initialized = false;

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
	hasReference = false;

	//soft clipped data
	softClippedEntry = 0;
	softClippedLength = NULL;
	softClippedBase = NULL;
	softClippedQuality = NULL;
}

TAlignment::TAlignment(unsigned int MaxSize){
	TAlignment();
	initStorage();

}

void TAlignment::clear(){
	empty = true;
	parsed = false;
	recalibrated = false;
	changed = false;
	passedFilters = false;
	quality = qualityOriginal;
}

void TAlignment::initStorage(){
	clear();

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
	softClippedEntry = 0;
	softClippedLength = new int[2];
	softClippedBase = new char*[2];
	softClippedQuality = new char*[2];
	softClippedBase[0] = new char[maxSize];
	softClippedBase[1] = new char[maxSize];
	softClippedQuality[0] = new char[maxSize];
	softClippedQuality[1] = new char[maxSize];

	initialized = true;

}

void TAlignment::fill(BamTools::BamAlignment & bamAlignment, int ReadGroupId){
	//clear
	clear();

	//add basic info
	chrNumber = bamAlignment.RefID;
	position = bamAlignment.Position;
	isReverseStrand = bamAlignment.IsReverseStrand();
	isProperPair = bamAlignment.IsProperPair();
	readGroupId = ReadGroupId;
}

void TAlignment::fillReferenceSequence(TFastaBuffer* fastaBuffer){
	if(!hasReference) //is this check really necessary?
		throw "No reference provided!";

	fastaBuffer->fill(chrNumber, position, position + alignedPos[length-1], referenceSequence);
};

void TAlignment::setFiltersPassed(bool passed){
	passedFilters = passed;
};

void TAlignment::setDistancesFromEnds(){
	//is it paired-end?
	if(isProperPair){
		if(isReverseStrand){
			//reverse (can be either first or second mate, but it's the one that comes second in bam file)
			//hence distance from 3' is given by f(dist since beginning of fragment) = f(insert - len + pos)
			//and distance from 5' is given as f(end of fragment) = f(len - pos - 1)
			int k = abs(bamAlignment.InsertSize) - length;
			int p = length - 1;
			for(int d=0; d<length; ++d){
				distFrom5Prime[d] = p - d;
				distFrom3Prime[d] = k + d;
			}
		} else {
			//forward (can be either first or second mate, but it's the one that comes first in bam file)
			//Hence distance from 3' is given as a function of pos
			//And distance from 5' is given by (length of fragment) - pos -1
			//NOTE! we ignore indels when calculating distance from 5' since we can not know this info.
			//Luckily, this has only minimal effect since these distances are far from fragment ends
			int p = abs(bamAlignment.InsertSize) - 1;
			for(int d=0; d<length; ++d){
				distFrom5Prime[d] = d;
				distFrom3Prime[d] = p - d;
			}
		}
	} else {
		//treat as single end
		int p = length - 1;
		if(isReverseStrand){
			//not in pair & reverse
			//Hence distance from 3' is just pos
			//And distance from 5' is just len - pos - 1
			for(int d=0; d<length; ++d){
				distFrom5Prime[d] = p - d;
				distFrom3Prime[d] = d;
			}

		} else {
			//not in pair & forward
			//Hence distance from 5' is just pos
			//And distance from 3' is given by len - pos - 1
			for(int d=0; d<length; ++d){
				distFrom5Prime[d] = d;
				distFrom3Prime[d] = p - d;
			}
		}
	}
};

void TAlignment::parse(bool & applyQualityFilter, bool & trimReads, int & minQual, int & maxQual, int & trimmingLength3Prime, int & trimmingLength5Prime, TGenotypeMap & genoMap, TQualityMap & qualityMap){
	if(!parsed){

		//first parse bases and qualities
		parseBasesQualities(genoMap, qualityMap);

		//then update distances from ends
		setDistancesFromEnds();

		//fill context for each base
		fillContext(genoMap);

		//apply filters
		//TODO: should that be on the recalibrated quality scores instead???
		if(applyQualityFilter)
			filterForBaseQuality(minQual, maxQual);
		if(trimReads)
			trimRead(trimmingLength3Prime, trimmingLength5Prime);

		parsed = true;
	}
};


void TAlignment::parseBasesQualities(TGenotypeMap & genoMap, TQualityMap & qualityMap){
	// iterate over CigarOps
	int d = 0; //index regarding data structures
	int k = 0; //index inside read
	int p = 0; //index regarding reference position (!= k for indels)
	softClippedEntry = 0; //softclipped bases to be added
	softClippedLength[0] = 0; softClippedLength[1] = 0;

	std::vector<BamTools::CigarOp>::const_iterator cigarIter = bamAlignment.CigarData.begin();
	std::vector<BamTools::CigarOp>::const_iterator cigarEnd  = bamAlignment.CigarData.end();
	int counter = 0;
	for( ; cigarIter != cigarEnd; ++cigarIter, ++counter ){
		const BamTools::CigarOp& op = (*cigarIter);
		switch ( op.Type ) {

			// for 'M', '=' or 'X': just copy
			case (BamTools::Constants::BAM_CIGAR_MATCH_CHAR)    :
			case (BamTools::Constants::BAM_CIGAR_SEQMATCH_CHAR) :
			case (BamTools::Constants::BAM_CIGAR_MISMATCH_CHAR) :
				for(unsigned int i=0; i<op.Length; ++i, ++d, ++k, ++p){
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
				for(unsigned int i=0; i<op.Length; ++i, ++d, ++k, ++p){
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
				for(unsigned int i=0; i<op.Length; ++i, ++d, ++k){
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
				for(unsigned int i=0; i<op.Length; ++i, ++d, ++k, ++p){
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
		throw "The lengths of the alignment and the quality scores of read '" + bamAlignment.Name + "' do not match!";
};



void TAlignment::fillContext(TGenotypeMap & genoMap){
	if(isReverseStrand){
		//reverse
		for(int d=0; d<(length-1); ++d){
//			std::cout << "getting Context for " << base[d-1] << ", " << base[d] << std::flush;
			//std::cout << " -> " << genoMap.contextMap[base[d-1]][base[d]] << std::endl;
			context[d] = genoMap.contextMap[base[d+1]][base[d]];
		}
		context[length-1] = genoMap.contextMap[N][base[length-1]];
	} else {
		//forward
		context[0] = genoMap.contextMap[N][base[0]];
		for(int d=1; d<length; ++d){
//			std::cout << "getting Context for " << base[d-1] << ", " << base[d] << std::flush;
			//std::cout << " -> " << genoMap.contextMap[base[d-1]][base[d]] << std::endl;
			context[d] = genoMap.contextMap[base[d-1]][base[d]];
		}
	}
};

void TAlignment::fillPmdProbabilities(TPMD* pmdObjects){
	for(int d=0; d<length; ++d){
		pmdCT[d] = pmdObjects[readGroupId].getProbCT(distFrom5Prime[d]);
		pmdGA[d] = pmdObjects[readGroupId].getProbGA(distFrom3Prime[d]);
	}
};

//--------------------------------------
//functions to access and modify data
//--------------------------------------
void TAlignment::filterForBaseQuality(int & minQual, int & maxQual){
	//set base to N if outside quality filter
	for(int d=0; d<length; ++d){
		if(qualityOriginal[d] < minQual || qualityOriginal[d] > maxQual){
			base[d] = N;
		}
	}
};

void TAlignment::filterForPrintingBaseQuality(std::string & qual, int & minQualForPrinting, int & maxQualForPrinting){
	//set base to N if outside quality filter
	for(std::string::iterator stringIt = qual.begin() ; stringIt < qual.end(); ++stringIt){
		if((int) *stringIt < minQualForPrinting)
			*stringIt = (char) minQualForPrinting;
		else if((int) *stringIt > maxQualForPrinting)
			*stringIt = (char) maxQualForPrinting;
	}
};

void TAlignment::trimRead(int & trimmingLength3Prime, int & trimmingLength5Prime){
	//set base to N at ends of read
	if(isReverseStrand){
		//distance from 3' is just pos
		//distance from 5' is len - pos - 1
		for(int d=0; d<trimmingLength3Prime; ++d)
			base[d] = N;
		for(int d=0; d<trimmingLength5Prime; ++d)
			base[length - 1 - d] = N;
	} else {
		//distance from 3' is len - pos - 1
		//distance from 5' is just pos
		for(int d=0; d<trimmingLength3Prime; ++d)
			base[length - 1 - d] = N;
		for(int d=0; d<trimmingLength5Prime; ++d)
			base[d] = N;
	}
};
void TAlignment::recalibrate(TRecalibration & recalObject, TQualityMap & qualityMap){
	//make sure read is parsed
	parse();

	if(recalObject.recalibrationChangesQualities()){
		//recalibrate quality scores
		for(int d=0; d<length; ++d){
			int k = length - d - 1;
			errorRates[d] = recalObject.getErrorRate(readGroupId, qualityOriginal[d], d, k, context[d]);
			qualityRecalibrated[d] = qualityMap.errorToQuality(errorRates[d]);
		}

		quality = qualityRecalibrated;
		changed = true;
	} else changed = false;
	recalibrated = true;

};

void TAlignment::recalibrate(TRecalibration & recalObject, TPMD* pmdObjects, TFastaBuffer* fastaBuffer, TQualityMap & qualityMap){
	//make sure read is parsed
	parse();

	//get reference sequence
	fillReferenceSequence(fastaBuffer);

	//get PMD probs
	fillPmdProbabilities(pmdObjects);

	//recalibrate quality scores
	for(int d=0; d<length; ++d){
		int k = length - d - 1;
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

void TAlignment::binQualityScores(TQualityMap & qualityMap){
	//make sure read is parsed
	parse();

	//bin quality scores as done by Illumina
	for(int d=0; d<length; ++d){
		quality[d] = qualityMap.illuminaQualityBins[quality[d]];
	}
	changed = true;
};

void TAlignment::addToPMDTables(TPMDTables & pmdTables, TFastaBuffer* fastaBuffer, TGenotypeMap & genoMap){
	//make sure read is parsed
	parse();

	//get reference sequence
	fillReferenceSequence(fastaBuffer);

	//tmp variables
	Base ref, read;

	//check if it is forward or reverse strand!
	if(isReverseStrand){
		for(int d=0; d<length; ++d){
			if(aligned[d] && base[d] != N){
				ref = genoMap.flipBase(referenceSequence[alignedPos[d]]);
				read = genoMap.baseToFlippedBase[base[d]];
				pmdTables.addForward(readGroupId, distFrom3Prime[d], ref, read);
				pmdTables.addReverse(readGroupId, distFrom5Prime[d], ref, read);
			}
		}
	} else {
		for(int d=0; d<length; ++d){
			if(aligned[d] && base[d] != N){
				ref = genoMap.getBase(referenceSequence[alignedPos[d]]);
				pmdTables.addForward(readGroupId, distFrom5Prime[d], ref, base[d]);
				pmdTables.addReverse(readGroupId, distFrom3Prime[d], ref, base[d]);
			}
		}
	}
};


double TAlignment::calculatePMDS(double & pi, TPMD* pmdObjects, TFastaBuffer* fastaBuffer){
	//make sure read is parsed
	parse();

	//variables
	double PMDS = 0.0;
	double probPMD, probNoPMD;
	double epsThird;
	double fourEpsThird;

	//get reference
	fillReferenceSequence(fastaBuffer);

	//get PMD probs
	fillPmdProbabilities(pmdObjects);

	//go over all bases in read
	for(int d=0; d<length; ++d){
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

void TAlignment::assessSoftClipping(int & S_left, int & middle, int & S_right){
	//count S, not S, S pattern from cigar string
	S_left = 0;
	S_right = 0;
	middle = 0;
	static bool reachedMiddle = false;

	std::vector<BamTools::CigarOp>::const_iterator cigarIter = bamAlignment.CigarData.begin();
	std::vector<BamTools::CigarOp>::const_iterator cigarEnd  = bamAlignment.CigarData.end();

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

void TAlignment::addToQualityTable(TQualityTable & qualTable){
	parse();
	qualTable.add(quality, length);
};



//--------------------------------------------
//functions to modify alignment
//--------------------------------------------
void TAlignment::updateOptionalSamField(std::string tag, float value){
	if(bamAlignment.HasTag(tag) == false) bamAlignment.AddTag(tag, "f", value);
	else bamAlignment.EditTag(tag, "f", value);
};

void TAlignment::downsampleAlignment(double& fraction, TRandomGenerator& randomGenerator){
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
void TAlignment::save(BamTools::BamWriter & bamWriter, TGenotypeMap & genoMap, int & minQualForPrinting, int & maxQualForPrinting){
	if(changed){
		//means that read has been modified.
		//Currently quality recalibration is the only possible change.
		//But will need to think how to deal with merging and such...

		//assume that only bases and quality scores where changed
		std::string tmpString, tmpString2;
//		tmpString.clear();
//		tmpString2.clear();
		for(int d=0; d<length; ++d){
			tmpString += genoMap.baseToChar[base[d]];
			tmpString2 += (char) quality[d];
		}
		bamAlignment.QueryBases = tmpString;
		bamAlignment.Qualities = tmpString2;
	}

	//make sure quality are printed within range
	filterForPrintingBaseQuality(bamAlignment.Qualities, minQualForPrinting, maxQualForPrinting);

	//now write alignment
	if(!bamWriter.SaveAlignment(bamAlignment))
		throw "Read '" + bamAlignment.Name + "' could not be written!";
};

void TAlignment::print(TGenotypeMap & genoMap){
	std::cout << "NAME:\t" << bamAlignment.Name << std::endl;
	std::cout << "LEN:\t" << length << std::endl;

	//print bases
	std::cout << "SEQ:\t";
	for(int d=0; d<length; ++d)
		std::cout << genoMap.getBaseAsChar(base[d]);
	std::cout << std::endl;

	//print qualities
	std::cout << "QUAL:\t";
	for(int d=0; d<length; ++d)
		std::cout << (char) quality[d];
	std::cout << std::endl;

	//print aligned pos
	std::cout << "POS:\t";
	for(int d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		if(aligned[d])
			std::cout << alignedPos[d];
		else
			std::cout << "-";
	}
	std::cout << std::endl;

	//print dist from 3'
	std::cout << "dist 3':\t";
	for(int d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		std::cout << distFrom3Prime[d];
	}
	std::cout << std::endl;

	//print dist from 5'
	std::cout << "dist 5':\t";
	for(int d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		std::cout << distFrom5Prime[d];
	}
	std::cout << std::endl;
}
