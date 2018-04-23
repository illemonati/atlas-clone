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
	mappingQuality = 0;
	passedFilters = false;
	parsed = false;
	changed = false;
	storageInitialized = false;
	recalibrated = false;
	hasReference = false;

	bases = NULL;
	aligned = NULL;
	alignedPos = NULL;
	qualityRecalibrated = NULL;
	quality = NULL;
	qualityOriginal = NULL;
	baseAsChar = NULL;
	softClippedEntry = 0;
	softClippedLength = NULL;
	softClippedBase = NULL;
	softClippedQuality = NULL;

	//per base data
/*
	base = NULL;
	context = NULL;
	errorRates = NULL;

	distFrom3Prime = NULL;
	distFrom5Prime = NULL;
	pmdCT = NULL;
	pmdGA = NULL;

	//soft clipped data

*/
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
	bases = new TBase[maxSize];

/*
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
*/

	storageInitialized = true;

}

void TAlignment::freeStorage(){
	if(storageInitialized){
		delete[] bases;
/*
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

		delete[] softClippedLength;
		delete[] softClippedBase[0];
		delete[] softClippedBase[1];
		delete[] softClippedBase;
		delete[] softClippedQuality[0];
		delete[] softClippedQuality[1];
		delete[] softClippedQuality;
*/

	}
	storageInitialized = false;
}

void TAlignment::fill(BamTools::BamAlignment & bamAlignment, int ReadGroupId){
	//clear
	clear();

	//add basic info
	chrNumber = bamAlignment.RefID;
	position = bamAlignment.Position;
	isReverseStrand = bamAlignment.IsReverseStrand();
	isProperPair = bamAlignment.IsProperPair();
	mappingQuality = bamAlignment.MapQuality;
	readGroupId = ReadGroupId;
}

void TAlignment::setFiltersPassed(bool passed){
	passedFilters = passed;
};

void TAlignment::setReferenceAdded(){
	hasReference = true;
}

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
				bases[d].posInRead = p - d;
				bases[d].posInReadRev = k + d;
			}
		} else {
			//forward (can be either first or second mate, but it's the one that comes first in bam file)
			//Hence distance from 3' is given as a function of pos
			//And distance from 5' is given by (length of fragment) - pos -1
			//NOTE! we ignore indels when calculating distance from 5' since we can not know this info.
			//Luckily, this has only minimal effect since these distances are far from fragment ends
			int p = abs(bamAlignment.InsertSize) - 1;
			for(int d=0; d<length; ++d){
				bases[d].posInRead = d;
				bases[d].posInReadRev = p - d;
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
				bases[d].posInRead = p - d;
				bases[d].posInReadRev = d;
			}

		} else {
			//not in pair & forward
			//Hence distance from 5' is just pos
			//And distance from 3' is given by len - pos - 1
			for(int d=0; d<length; ++d){
				bases[d].posInRead = d;
				bases[d].posInReadRev = p - d;
			}
		}
	}
};

void TAlignment::parse(TGenotypeMap & genoMap, TQualityMap & qualityMap){
	if(!parsed){

		if(!storageInitialized) throw "Alignment storage was not initialized!";

		//first parse bases and qualities
		parseBasesQualities(genoMap, qualityMap);

		//then update distances from ends
		setDistancesFromEnds();

		//fill context for each base
		fillContext(genoMap);

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
					bases[d].base = genoMap.getBase(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					qualityOriginal[d] = (int) bamAlignment.Qualities[k];
					bases[d].errorRate = qualityMap.qualityToErrorMap[(int) bamAlignment.Qualities[k]];
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
					bases[d].base = genoMap.getBase(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					qualityOriginal[d] = -1;
					bases[d].errorRate = qualityMap.qualityToErrorMap[(int) bamAlignment.Qualities[k]];
					aligned[d] = false;
					alignedPos[d] = -1;
				}
				break;

			//for 'I' - insertion: copy bases, but put aligned pos to
			case (BamTools::Constants::BAM_CIGAR_INS_CHAR)      :
				for(unsigned int i=0; i<op.Length; ++i, ++d, ++k){
					bases[d].base = genoMap.getBase(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					qualityOriginal[d] = (int) (char) bamAlignment.Qualities[k];
					bases[d].errorRate = qualityMap.qualityToErrorMap[(int) bamAlignment.Qualities[k]];
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
					bases[d].base = genoMap.getBase(bamAlignment.QueryBases[k]);
					baseAsChar[d] = bamAlignment.QueryBases[k];
					qualityOriginal[d] = (int) bamAlignment.Qualities[k];
					bases[d].errorRate = qualityMap.qualityToErrorMap[(int) bamAlignment.Qualities[k]];
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
//			std::cout << "getting Context for " << base[d-1] << ", " << bases[d].base << std::flush;
			//std::cout << " -> " << genoMap.contextMap[base[d-1]][bases[d].base] << std::endl;
			bases[d].context = genoMap.contextMap[bases[d+1].base][bases[d].base];
		}
		bases[length-1].context = genoMap.contextMap[N][bases[length-1].base];
	} else {
		//forward
		bases[0].context = genoMap.contextMap[N][bases[0].base];
		for(int d=1; d<length; ++d){
//			std::cout << "getting Context for " << base[d-1] << ", " << bases[d].base << std::flush;
			//std::cout << " -> " << genoMap.contextMap[base[d-1]][bases[d].base] << std::endl;
			bases[d].context = genoMap.contextMap[bases[d-1].context][bases[d].context];
		}
	}
};

void TAlignment::fillPmdProbabilities(TPMD* pmdObjects){
	for(int d=0; d<length; ++d){
		bases[d].PMD_CT = pmdObjects[readGroupId].getProbCT(bases[d].posInRead);
//		bases[d].PMD_CT = pmdObjects[readGroupId].getProbCT(distFrom5Prime[d]);
//		bases[d].PMD_GA = pmdObjects[readGroupId].getProbGA(distFrom3Prime[d]);
	}
};

//--------------------------------------
//functions to access and modify data
//--------------------------------------
void TAlignment::filterForBaseQuality(int & minQual, int & maxQual){
	//set base to N if outside quality filter
	for(int d=0; d<length; ++d){
		if(qualityOriginal[d] < minQual || qualityOriginal[d] > maxQual){
			bases[d].base = N;
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
			bases[d].base = N;
		for(int d=0; d<trimmingLength5Prime; ++d)
			bases[length - 1 - d].base = N;
	} else {
		//distance from 3' is len - pos - 1
		//distance from 5' is just pos
		for(int d=0; d<trimmingLength3Prime; ++d)
			bases[length - 1 - d].base = N;
		for(int d=0; d<trimmingLength5Prime; ++d)
			bases[d].base = N;
	}
};
void TAlignment::recalibrate(TRecalibration & recalObject, TQualityMap & qualityMap){
	//make sure read is parsed and has reference
	if(!parsed) throw "Read was not parsed!";

	if(recalObject.recalibrationChangesQualities()){
		//recalibrate quality scores
		for(int d=0; d<length; ++d){
			int k = length - d - 1;
			bases[d].errorRate = recalObject.getErrorRate(readGroupId, qualityOriginal[d], d, k, bases[d].context);
			qualityRecalibrated[d] = qualityMap.errorToQuality(bases[d].errorRate);

//			bases[d].error = recalObject.getErrorRate(readGroupId, qualityOriginal[d], d, k, bases[d].context);
//			qualityRecalibrated[d] = qualityMap.errorToQuality(bases[d].error);
		}

		quality = qualityRecalibrated;
		changed = true;
	} else changed = false;
	recalibrated = true;

};

void TAlignment::recalibrate(TRecalibration & recalObject, TPMD* pmdObjects, TFastaBuffer* fastaBuffer, TQualityMap & qualityMap){
	//make sure read is parsed and has reference
	if(!parsed) throw "Read was not parsed!";
	if(!hasReference) throw "Reference was not added!";

	//get PMD probs
	fillPmdProbabilities(pmdObjects);

	//recalibrate quality scores
	for(int d=0; d<length; ++d){
		int k = length - d - 1;
		if(recalObject.recalibrationChangesQualities())
			bases[d].errorRate = recalObject.getErrorRate(readGroupId, qualityOriginal[d], d, k, bases[d].context);

		//now add effect of PMD
		if(aligned[d]){
			if(bases[d].base == T && referenceSequence[alignedPos[d]] == 'C')
				bases[d].errorRate = 1.0 - ((1.0 - bases[d].errorRate)*(1.0 - bases[d].PMD_CT)); //this is mapDamage2, Krishna: qual*(1-pmdCT) + (1-qual)*pmdCT;
			else if(bases[d].base == A && referenceSequence[alignedPos[d]] == 'G')
				bases[d].errorRate = 1.0 - ((1.0 - bases[d].errorRate)*(1.0 - bases[d].PMD_GA)); //this is mapDamage2, Krishna: qual*(1-pmdGA) + (1-qual)*pmdGA;
		}

		qualityRecalibrated[d] = qualityMap.errorToQuality(bases[d].errorRate);
	}

	//set pointer to recalibrated scores
	quality = qualityRecalibrated;
	recalibrated = true;
	changed = true;
};

void TAlignment::binQualityScores(TQualityMap & qualityMap){
	//make sure read is parsed
	if(!parsed) throw "Read was not parsed!";

	//bin quality scores as done by Illumina
	for(int d=0; d<length; ++d){
		quality[d] = qualityMap.illuminaQualityBins[quality[d]];
	}
	changed = true;
};

void TAlignment::addToPMDTables(TPMDTables & pmdTables, TGenotypeMap & genoMap){
	//make sure read is parsed and has reference
	if(!parsed) throw "Read was not parsed!";
	if(!hasReference) throw "Reference was not added!";

	//tmp variables
	Base ref, read;

	//check if it is forward or reverse strand!
	if(isReverseStrand){
		for(int d=0; d<length; ++d){
			if(aligned[d] && bases[d].base != N){
				ref = genoMap.flipBase(referenceSequence[alignedPos[d]]);
				read = genoMap.baseToFlippedBase[bases[d].base];
				pmdTables.addForward(readGroupId, bases[d].posInReadRev, ref, read);
				pmdTables.addReverse(readGroupId, bases[d].posInRead, ref, read);
			}
		}
	} else {
		for(int d=0; d<length; ++d){
			if(aligned[d] && bases[d].base != N){
				ref = genoMap.getBase(referenceSequence[alignedPos[d]]);
				pmdTables.addForward(readGroupId, bases[d].posInRead, ref, bases[d].base);
				pmdTables.addReverse(readGroupId, bases[d].posInReadRev, ref, bases[d].base);
			}
		}
	}
};


double TAlignment::calculatePMDS(double & pi, TPMD* pmdObjects){
	//make sure read is parsed and has reference
	if(!parsed) throw "Read was not parsed!";
	if(!hasReference) throw "Reference was not added!";

	//variables
	double PMDS = 0.0;
	double probPMD, probNoPMD;
	double epsThird;
	double fourEpsThird;

	//get PMD probs
	fillPmdProbabilities(pmdObjects);

	//go over all bases in read
	for(int d=0; d<length; ++d){
		//limit to aligned positions
		if(aligned[d]){
			//Prepare variables
			epsThird = bases[d].errorRate / 3.0;
			fourEpsThird = 4.0 * epsThird;

			//calc likelihoods
			if(referenceSequence[alignedPos[d]] == 'A'){
				if(bases[d].base == A){
					probPMD = 1.0 - bases[d].errorRate - pi + fourEpsThird*pi + bases[d].PMD_GA*pi/3.0*(1.0-fourEpsThird);
					probNoPMD = 1.0 - bases[d].errorRate - pi + fourEpsThird*pi;
				}
				else if(bases[d].base == C){ //ok
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi - pi*bases[d].PMD_CT*(fourEpsThird-1.0);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else if(bases[d].base == G){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + pi*bases[d].PMD_GA*(fourEpsThird-1.0);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else if(bases[d].base == T){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + pi*bases[d].PMD_CT*(1.0-fourEpsThird);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
			}
			else if (referenceSequence[alignedPos[d]] == 'C'){
				if(bases[d].base == A){
					probPMD = bases[d].errorRate + pi - 2.0*bases[d].errorRate*pi + pi*bases[d].PMD_GA*(1.0-fourEpsThird);
					probNoPMD = bases[d].errorRate + pi - 2.0*bases[d].errorRate*pi;
				}
				else if(bases[d].base == C){
					probPMD = 1.0 - pi - bases[d].errorRate + fourEpsThird*pi + (1.0-pi)*bases[d].PMD_CT*(fourEpsThird-1.0);
					probNoPMD = 1.0 - pi - bases[d].errorRate + fourEpsThird*pi;
				}
				else if(bases[d].base == G){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + bases[d].PMD_GA*(fourEpsThird*pi - pi);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else if(bases[d].base == T){
					probPMD = epsThird + (1.0-pi)*bases[d].PMD_CT*(1.0-fourEpsThird);
					probNoPMD = epsThird;
				}
			}
			else if (referenceSequence[alignedPos[d]] == 'G'){
				if(bases[d].base == A){
					probPMD = bases[d].PMD_GA*(3.0-3.0*pi+4.0*bases[d].errorRate+4.0*bases[d].errorRate*pi) + bases[d].errorRate - fourEpsThird*pi + pi;
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else if(bases[d].base == C){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + pi*bases[d].PMD_CT*(fourEpsThird - 1.0);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else if(bases[d].base == G){
					probPMD = 1.0 - bases[d].errorRate - pi + fourEpsThird*pi + (1.0-pi)*bases[d].PMD_GA*(fourEpsThird-1.0);
					probNoPMD = 1.0 - bases[d].errorRate - pi + fourEpsThird*pi;
				}
				else if(bases[d].base == T){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + pi*bases[d].PMD_CT*(1.0-fourEpsThird);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
			}
			else if(referenceSequence[alignedPos[d]] == 'T'){
				if(bases[d].base == A){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi - epsThird*pi*bases[d].PMD_CT + pi*bases[d].PMD_GA*(1.0-bases[d].errorRate);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else if(bases[d].base == C){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + pi*bases[d].PMD_CT*(fourEpsThird - 1.0);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else if(bases[d].base == G){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + pi*bases[d].PMD_GA*(fourEpsThird - 1.0);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else if(bases[d].base == T){
					probPMD = 1.0 - bases[d].errorRate - pi + fourEpsThird*pi + bases[d].PMD_CT*(pi/3.0-fourEpsThird*pi/3.0);
					probNoPMD = 1.0 - bases[d].errorRate - pi + fourEpsThird*pi;
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
	//make sure read is parsed
	if(!parsed) throw "Read was not parsed!";
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
			bases[d].base = N;
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
			tmpString += genoMap.baseToChar[bases[d].base];
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
		std::cout << genoMap.getBaseAsChar(bases[d].base);
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
		std::cout << bases[d].posInReadRev;
	}
	std::cout << std::endl;

	//print dist from 5'
	std::cout << "dist 5':\t";
	for(int d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		std::cout << bases[d].posInRead;
	}
	std::cout << std::endl;
}
