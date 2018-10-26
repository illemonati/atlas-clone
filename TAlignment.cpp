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
	lastAlignedPos = 0;
	lastPositionPlusOne = 0;
	isReverseStrand = false;
	isProperPair = false;
	mappingQuality = 0;
	passedFilters = false;
	parsed = false;
	changed = false;
	storageInitialized = false;
	recalibrated = false;
	std::string referenceSequence = "";
	hasReference = false;

	bases = NULL;
//	quality = NULL;
	qualityOriginal = NULL;
	numInsertions = -1;
	numDeletions = -1;
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
	maxSize = MaxSize;
	initStorage();
}

TAlignment::TAlignment(TAlignment & Alignment){
	//details
	empty = true;
	maxSize = Alignment.maxSize;
	alignmentName = Alignment.alignmentName;
	length = Alignment.length;
	chrNumber = Alignment.chrNumber;
	readGroupId = Alignment.readGroupId;
	readGroup = Alignment.readGroup;
	position = Alignment.position;
	lastAlignedPos = Alignment.lastAlignedPos;
	lastPositionPlusOne = Alignment.lastPositionPlusOne;
	isReverseStrand = Alignment.isReverseStrand;
	isProperPair = Alignment.isProperPair;
	mappingQuality = Alignment.mappingQuality;
	passedFilters = Alignment.passedFilters;
	parsed = Alignment.parsed;
	changed = Alignment.changed;
	storageInitialized = Alignment.storageInitialized;
	recalibrated = Alignment.recalibrated;
	hasReference = Alignment.hasReference;
	referenceSequence = Alignment.referenceSequence;
	qualityOriginal = Alignment.qualityOriginal;
	numInsertions = Alignment.numInsertions;
	numDeletions = Alignment.numDeletions;
	softClippedEntry = Alignment.softClippedEntry;
	softClippedLength = Alignment.softClippedLength;
	softClippedBase = Alignment.softClippedBase;
	softClippedQuality = Alignment.softClippedQuality;
	bases = Alignment.bases;
}

void TAlignment::clear(){
	position = -1;
	alignmentName = "empty";
	empty = true;
	parsed = false;
	recalibrated = false;
	changed = false;
	passedFilters = false;
//	quality = qualityOriginal;
}

void TAlignment::initStorage(){
	clear();
	bases = new TBase[maxSize];

	//other
	qualityOriginal = new int[maxSize];

	//soft clipped data
	softClippedLength = new int[2];
	softClippedBase = new char*[2];
	softClippedQuality = new char*[2];
	softClippedBase[0] = new char[maxSize];
	softClippedBase[1] = new char[maxSize];
	softClippedQuality[0] = new char[maxSize];
	softClippedQuality[1] = new char[maxSize];

	storageInitialized = true;

}

void TAlignment::freeStorage(){
	if(storageInitialized){
		delete[] bases;
		delete[] qualityOriginal;

		delete[] softClippedLength;
		delete[] softClippedBase[0];
		delete[] softClippedBase[1];
		delete[] softClippedBase;
		delete[] softClippedQuality[0];
		delete[] softClippedQuality[1];
		delete[] softClippedQuality;

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


*/

	}
	storageInitialized = false;
}

void TAlignment::fill(BamTools::BamAlignment & BamAlignment, int ReadGroupId){
	//clear
	clear();

	//add basic info
	bamAlignment = BamAlignment;
	chrNumber = bamAlignment.RefID;
	position = bamAlignment.Position;
	alignmentName = bamAlignment.Name;
	isReverseStrand = bamAlignment.IsReverseStrand();
	isProperPair = bamAlignment.IsProperPair();
	mappingQuality = bamAlignment.MapQuality;
	readGroupId = ReadGroupId;

	empty = false;
}

void TAlignment::setReferenceAdded(){
	hasReference = true;
}

void TAlignment::setDistancesFromEnds(){
	//is it paired-end?
	if(isProperPair){
		if(isReverseStrand){
			//reverse (can be either first or second mate, but it's the one that comes second in bam file)
			//hence distance from 5' is given by f(dist since beginning of fragment) = f(insert - len + pos)
			//and distance from 3' is given as f(end of fragment) = f(len - pos - 1)
			int k = abs(bamAlignment.InsertSize) - (length - softClippedLength[1]) + numInsertions - numDeletions;
			int p = length - 1 - softClippedLength[1];
			for(int d=0; d<length; ++d){
				bases[d].posInRead = p - d; //dist from 5'
				bases[d].posInReadRev = k + d; //dist from 3'
			}
		} else {
			//forward (can be either first or second mate, but it's the one that comes first in bam file)
			//Hence distance from 5' is given as a function of pos
			//And distance from 3' is given by (length of fragment) - pos -1
			//NOTE! we ignore indels when calculating distance from 5' since we can not know this info.
			//Luckily, this has only minimal effect since these distances are far from fragment ends
			int p = abs(bamAlignment.InsertSize) - 1 + numInsertions - numDeletions;
			for(int d=0; d<length; ++d){
				bases[d].posInRead = d - softClippedLength[0]; //dist from 5'
				bases[d].posInReadRev = p - d + softClippedLength[0]; //dist from 3'
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
				bases[d].posInRead = p - d - softClippedLength[1]; //dist from 5'
				bases[d].posInReadRev = d - softClippedLength[0]; //dist from 3'
			}

		} else {
			//not in pair & forward
			//Hence distance from 5' is just pos
			//And distance from 3' is given by len - pos - 1
			for(int d=0; d<length; ++d){
				bases[d].posInRead = d - softClippedLength[0]; //dist from 5'
				bases[d].posInReadRev = p - d - softClippedLength[1]; //dist from 3'
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
	numInsertions = 0;
	numDeletions = 0;

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

				//soft-clipped bases on 5' are before bamAlignment.Position
				for(unsigned int i=0; i<op.Length; ++i, ++d, ++k, ++p){
					bases[d].base = genoMap.getBase(bamAlignment.QueryBases[k]);
					qualityOriginal[d] = (int) bamAlignment.Qualities[k];
					bases[d].errorRate = qualityMap.qualityToErrorMap[(int) bamAlignment.Qualities[k]];
					if(bases[d].base == N)
						bases[d].aligned = false;
					else {
						bases[d].aligned = true;
						bases[d].alignedPos = p;
					}
				}
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
				break;

			//for 'S' - soft clip: ignore bases, but increase k
			case (BamTools::Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
				//add bases to softclipped entries
				for(unsigned int i=0; i<op.Length; ++i, ++d, ++k){
					//soft-clipped bases on 5' are before bamAlignment.Position
					softClippedBase[softClippedEntry][softClippedLength[softClippedEntry]] = bamAlignment.QueryBases[k];
					softClippedQuality[softClippedEntry][softClippedLength[softClippedEntry]] = bamAlignment.Qualities[k];
					++softClippedLength[softClippedEntry];
					//need to initialize quality for quality filter and bases for context
					bases[d].base = genoMap.getBase(bamAlignment.QueryBases[k]);
					qualityOriginal[d] = -1;
					bases[d].errorRate = qualityMap.qualityToErrorMap[(int) bamAlignment.Qualities[k]];
					bases[d].aligned = false;
					bases[d].alignedPos = -1;
				}
				break;

			//for 'I' - insertion: copy bases, but put aligned pos to
			case (BamTools::Constants::BAM_CIGAR_INS_CHAR)      :
				for(unsigned int i=0; i<op.Length; ++i, ++d, ++k){
					bases[d].base = genoMap.getBase(bamAlignment.QueryBases[k]);
					qualityOriginal[d] = (int) (char) bamAlignment.Qualities[k];
					bases[d].errorRate = qualityMap.qualityToErrorMap[(int) bamAlignment.Qualities[k]];
					bases[d].aligned = false;
					++numInsertions;
					bases[d].alignedPos = -1;
				}
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
				break;


			// for 'D' - deletion: just add to postion
			case (BamTools::Constants::BAM_CIGAR_DEL_CHAR) :
				p += op.Length;
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
				++numDeletions;
				break;

			// for 'N' - skipped region: copy but say that bases were not aligned
			case (BamTools::Constants::BAM_CIGAR_REFSKIP_CHAR) :
				for(unsigned int i=0; i<op.Length; ++i, ++d, ++k, ++p){
					bases[d].base = genoMap.getBase(bamAlignment.QueryBases[k]);
					qualityOriginal[d] = (int) bamAlignment.Qualities[k];
					bases[d].errorRate = qualityMap.qualityToErrorMap[(int) bamAlignment.Qualities[k]];
					bases[d].aligned = false;
					bases[d].alignedPos = p;
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

	//update length and last aligned position
	length = k;
	lastPositionPlusOne = position + length;
	lastAlignedPos = p - 1;
	if(passedFilters && length != bamAlignment.Length)
		throw "The lengths of the alignment and the quality scores of read '" + bamAlignment.Name + "' do not match!";
};

void TAlignment::fillContext(TGenotypeMap & genoMap){
	if(isReverseStrand){
		//reverse
		for(int d=0; d<(length-1); ++d)
			bases[d].context = genoMap.contextMap[bases[d+1].base][bases[d].base];
		bases[length-1].context = genoMap.contextMap[N][bases[length-1].base];
	} else {
		//forward
		bases[0].context = genoMap.contextMap[N][bases[0].base];
		for(int d=1; d<length; ++d)
			bases[d].context = genoMap.contextMap[bases[d-1].base][bases[d].base];
	}
};

void TAlignment::fillPmdProbabilities(TPMD* pmdObjects){
	for(int d=0; d<length; ++d){
		bases[d].PMD_CT = pmdObjects[readGroupId].getProbCT(bases[d].posInRead);
		bases[d].PMD_GA = pmdObjects[readGroupId].getProbGA(bases[d].posInReadRev);
	}
};

//--------------------------------------
//functions to modify data
//--------------------------------------
void TAlignment::filterForBaseQuality(int & minQual, int & maxQual){
	//set base to N if outside quality filter
	for(int d=0; d<length; ++d){
		if(qualityOriginal[d] < minQual || qualityOriginal[d] > maxQual){
			bases[d].base = N;
		}
	}
	changed = true;
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
	//TODO: check if these distances correspond to those calculated in setDistances function
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
	changed = true;
};


/*
void TAlignment::recalibrate(TRecalibration & recalObject, TQualityMap & qualityMap){
	//make sure read is parsed and has reference
	if(!parsed) throw "Read was not parsed!";

	if(recalObject.recalibrationChangesQualities()){
		//recalibrate quality scores
		for(int d=0; d<length; ++d){
//			int k = length - d - 1;
			bases[d].errorRate = recalObject.getErrorRate(bases[d]);
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
//		int k = length - d - 1;
		if(recalObject.recalibrationChangesQualities())
			bases[d].errorRate = recalObject.getErrorRate(bases[d]);

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
*/


void TAlignment::fillReadGroupInfo(int & readGroupID){
	for(int d=0; d<length; ++d){
		bases[d].readGroup = readGroupID;
	}
}

void TAlignment::binQualityScores(TQualityMap & qualityMap){
	//make sure read is parsed
	if(!parsed) throw "Read was not parsed!";

	//bin quality scores as done by Illumina
	for(int d=0; d<length; ++d){
		bases[d].errorRate = qualityMap.qualityToError(qualityMap.illuminaQualityBins[qualityOriginal[d]]);
	}
	changed = true;
};

void TAlignment::updateOptionalSamField(std::string tag, float value){
	if(bamAlignment.HasTag(tag) == false) bamAlignment.AddTag(tag, "f", value);
	else bamAlignment.EditTag(tag, "f", value);
};

void TAlignment::updateOptionalSamField(std::string tag, std::string value){
	if(bamAlignment.HasTag(tag) == false) bamAlignment.AddTag(tag, "Z", value);
	else bamAlignment.EditTag(tag, "Z", value);
};

void TAlignment::downsampleAlignment(double& fraction, TRandomGenerator& randomGenerator, TQualityMap & qualMap){
	for(int d=0; d<length; ++d){
		double r = randomGenerator.getRand();
		if(r < fraction){
			bases[d].base = N;
			bases[d].errorRate = qualMap.qualityToError(0);
		}
	}
	changed = true;
}
//--------------------------------------
//functions to fill other classes
//--------------------------------------

void TAlignment::addToPMDTables(TPMDTables & pmdTables, TGenotypeMap & genoMap){
	//make sure read is parsed and has reference
	if(!parsed) throw "Read was not parsed!";
	if(!hasReference) throw "Reference was not added!";

	//tmp variables
	Base ref, read;

	//check if it is forward or reverse strand!
	if(isReverseStrand){
		for(int d=0; d<length; ++d){
			if(bases[d].aligned && bases[d].base != N){
				ref = genoMap.flipBase(referenceSequence[bases[d].alignedPos]);
				read = genoMap.baseToFlippedBase[bases[d].base];
				pmdTables.addForward(readGroupId, bases[d].posInReadRev, ref, read);
				pmdTables.addReverse(readGroupId, bases[d].posInRead, ref, read);
			}
		}
	} else {
		for(int d=0; d<length; ++d){
			if(bases[d].aligned && bases[d].base != N){
				ref = genoMap.getBase(referenceSequence[bases[d].alignedPos]);
				pmdTables.addForward(readGroupId, bases[d].posInRead, ref, bases[d].base);
				pmdTables.addReverse(readGroupId, bases[d].posInReadRev, ref, bases[d].base);
			}
		}
	}
};

void TAlignment::recalibrateWithPMD(TRecalibration* recalObject, TQualityMap & qualMap){
	//make sure read is parsed and has reference
	if(!parsed) throw "Read was not parsed!";
	if(!hasReference) throw "Reference was not added!";

	for(int d=0; d<length; ++d){
//		int k = length - d - 1;
		if(bases[d].aligned){
			//recalibrate quality scores
			if(recalObject->recalibrationChangesQualities())
				bases[d].errorRate = recalObject->getErrorRate(bases[d]);

			if(bases[d].context > 20) throw "there is a invalid context in alignment " + alignmentName + " at position " + toString(d);

			//now add effect of PMD
			if(bases[d].base == T && referenceSequence[bases[d].alignedPos] == 'C')
				bases[d].errorRate = 1.0 - ((1.0 - bases[d].errorRate)*(1.0 - bases[d].PMD_CT)); //this is mapDamage2, Krishna: qual*(1-pmdCT) + (1-qual)*pmdCT;
			else if(bases[d].base == A && referenceSequence[bases[d].alignedPos] == 'G')
				bases[d].errorRate = 1.0 - ((1.0 - bases[d].errorRate)*(1.0 - bases[d].PMD_GA)); //this is mapDamage2, Krishna: qual*(1-pmdGA) + (1-qual)*pmdGA;
		} else {
			bases[d].errorRate = qualMap.qualityToErrorMap[qualityOriginal[d]];
		}
	}

	recalibrated = true;
	changed = true;
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
		if(bases[d].aligned){
			//Prepare variables
			epsThird = bases[d].errorRate / 3.0;
			fourEpsThird = 4.0 * epsThird;

			//calc likelihoods
			if(referenceSequence[bases[d].alignedPos] == 'A'){
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
			else if (referenceSequence[bases[d].alignedPos] == 'C'){
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
			else if (referenceSequence[bases[d].alignedPos] == 'G'){
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
			else if(referenceSequence[bases[d].alignedPos] == 'T'){
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

int TAlignment::measureOverlap(){
	//make sure read is parsed
	if(!parsed) throw "Read was not parsed!";

	if(isProperPair && passedFilters){
		if(!isReverseStrand){
			int k = length - softClippedLength[1] - 1;
			while(bases[k].alignedPos < 0)
				--k;
			int endPos = position + bases[k].alignedPos;
			int overlap = endPos - bamAlignment.MatePosition;

			if(overlap < 0)
				//there is no overlap
				return 0;
			else
				return overlap;
		} else
			//not relevant
			return -1;

	} else
		//not relevant
		return -1;
}

void TAlignment::addToQualityTable(TQualityTable & qualTable, TQualityMap & qualMap){
	//make sure read is parsed
	if(!parsed) throw "Read was not parsed!";
	for(int d=0; d<length; ++d){
		const int qual = qualMap.errorToQuality(bases[d].errorRate);
		qualTable.add(qual);
	}
};





//--------------------------------------------
//functions to write / print alignment
//--------------------------------------------
void TAlignment::setToSingleEnd(){
	bamAlignment.Length = bamAlignment.AlignedBases.size();
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp(BamTools::Constants::BAM_CIGAR_MATCH_CHAR, bamAlignment.Length));
	bamAlignment.SetIsFirstMate(false);
	bamAlignment.SetIsSecondMate(false);
	bamAlignment.SetIsPaired(false);
	bamAlignment.SetIsProperPair(false);
	bamAlignment.SetIsMateReverseStrand(false);
	bamAlignment.SetIsReverseStrand(false); //the read that comes first in BAM is always fwd strand
	bamAlignment.MateRefID = -1;
	bamAlignment.MatePosition = -1;
	bamAlignment.InsertSize = 0;
}


void TAlignment::save(BamTools::BamWriter & bamWriter, TGenotypeMap & genoMap, int & minQualForPrinting, int & maxQualForPrinting, TQualityMap & qualMap){
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

			tmpString2 += (char) qualMap.errorToQuality(bases[d].errorRate);
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

void TAlignment::print(TGenotypeMap & genoMap, TQualityMap & qualMap){
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
		std::cout << (char) qualMap.errorToQuality(bases[d].errorRate);
	std::cout << std::endl;

	//print aligned pos
	std::cout << "POS:\t";
	for(int d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		if(bases[d].aligned)
			std::cout << bases[d].alignedPos;
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
