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
	fragmentLength = 0;
	chrNumber = 0;
	readGroupId = 0;
	position = 0;
	matePosition = 0;
	lastAlignedPos = 0;
	lastPositionPlusOne = 0;
	lastAlignedPositionWithRespectToRef = 0;
	insertSize = 0;
	isReverseStrand = false;
	isPaired = false;
	isProperPair = false;
	mappingQuality = 0;
	isSecondMate = false;
	parsed = false;
	changed = false;
	storageInitialized = false;
	recalibrated = false;
	std::string referenceSequence = "";
	hasReference = false;

	bases = NULL;
	numInsertions = -1;
	numDeletions = -1;
	softClippedEntry = 0;
	softClippedLength = NULL;
	softClippedBase = NULL;
	softClippedQuality = NULL;
};

TAlignment::TAlignment(unsigned int MaxSize){
	TAlignment();
	maxSize = MaxSize;
	initStorage();
};

TAlignment::TAlignment(const TAlignment & Alignment){
	TAlignment();

	//initialize storage
	maxSize = Alignment.maxSize;
	initStorage();
	storageInitialized = true;

	//copy content
	bamAlignment = Alignment.bamAlignment;
	alignmentName = Alignment.alignmentName;
	length = Alignment.length;
	chrNumber = Alignment.chrNumber;
	readGroupId = Alignment.readGroupId;
	position = Alignment.position;
	lastAlignedPos = Alignment.lastAlignedPos;
	lastPositionPlusOne = Alignment.lastPositionPlusOne;
	lastAlignedPositionWithRespectToRef = Alignment.lastAlignedPositionWithRespectToRef;
	isReverseStrand = Alignment.isReverseStrand;
	isPaired = Alignment.isPaired;
	isProperPair = Alignment.isProperPair;
	mappingQuality = Alignment.mappingQuality;
	parsed = Alignment.parsed;
	changed = Alignment.changed;
	recalibrated = Alignment.recalibrated;
	hasReference = Alignment.hasReference;
	referenceSequence = Alignment.referenceSequence;
	numInsertions = Alignment.numInsertions;
	numDeletions = Alignment.numDeletions;

	empty = false;

	//copy data from arrays
	softClippedEntry = Alignment.softClippedEntry;
	std::copy(Alignment.softClippedLength, Alignment.softClippedLength + 2, softClippedLength);
	std::copy(Alignment.softClippedBase[0], Alignment.softClippedBase[0] + Alignment.softClippedLength[0], softClippedBase[0]);
	std::copy(Alignment.softClippedBase[1], Alignment.softClippedBase[1] + Alignment.softClippedLength[1], softClippedBase[1]);

//	for(int e=0; e<2; ++e){
//		for(int b=0; b<Alignment.softClippedLength[e]; ++b){
//			std::cout << "Alignment.softClippedBase[e][b] " << std::flush;
//			std::cout << Alignment.softClippedBase[e][b] << std::endl;
//			softClippedBase[e][b] = Alignment.softClippedBase[e][b];
//			softClippedQuality[e][b] = Alignment.softClippedQuality[e][b];
//
//			std::cout << "this.softClippedBase[e][b] " << std::flush;
//			std::cout << Alignment.softClippedBase[e][b] << std::endl;
//		}
//	}
//	throw "done";

	std::copy(Alignment.bases, Alignment.bases + Alignment.maxSize, bases);
}

void TAlignment::clear(){
	position = -1;
	alignmentName = "empty";
	empty = true;
	parsed = false;
	recalibrated = false;
	changed = false;
};

void TAlignment::initStorage(){
	clear();

	//data
	bases = new TBase[maxSize];

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
		delete[] softClippedLength;
		delete[] softClippedBase[0];
		delete[] softClippedBase[1];
		delete[] softClippedBase;
		delete[] softClippedQuality[0];
		delete[] softClippedQuality[1];
		delete[] softClippedQuality;
	}
	storageInitialized = false;
};

void TAlignment::fill(BamTools::BamAlignment & BamAlignment, const uint16_t ReadGroupId){
	//clear
	clear();

	//add basic info
	bamAlignment = BamAlignment;
	chrNumber = bamAlignment.RefID;
	position = bamAlignment.Position;
	alignmentName = bamAlignment.Name;
	isReverseStrand = bamAlignment.IsReverseStrand();
	isPaired = bamAlignment.IsPaired();
	isProperPair = bamAlignment.IsProperPair();
	mappingQuality = bamAlignment.MapQuality;
	readGroupId = ReadGroupId;
	insertSize = bamAlignment.InsertSize;
	isSecondMate = bamAlignment.IsSecondMate();
	matePosition = bamAlignment.MatePosition;

	empty = false;
};

void TAlignment::setReferenceAdded(){
	hasReference = true;
}

void TAlignment::setDistancesFromEnds(){
	//Set distances in ORIGINAL FRAGMENT (i.e. 5' end is where sequencing started, NOT how it aligns to reference)
	//is it paired-end?
	if(isProperPair){
		int fragmentLen = abs(insertSize) + numInsertions - numDeletions;
		if(isReverseStrand){
			//reverse (can be either first or second mate, but it's the one that comes second in bam file)
			//and distance from 5' is given as f(end of fragment) = f(len - pos - 1)
			//hence distance from 3' is given by f(dist since beginning of fragment) = f(insert - len + pos)
			int k = abs(fragmentLen) - (length - softClippedLength[1]);
			int l = length - 1 - softClippedLength[1];
			for(int pos=0; pos<length; ++pos){
				bases[pos].data.distFrom5Prime = l - pos; //dist from 5'
				bases[pos].data.distFrom3Prime = k + pos; //dist from 3'
			}
		} else {
			//forward (can be either first or second mate, but it's the one that comes first in bam file)
			//Hence distance from 5' is given as a function of pos
			//And distance from 3' is given by (length of fragment) - pos -1
			//NOTE! we ignore indels when calculating distance from 5' since we can not know this info.
			//Luckily, this has only minimal effect since these distances are far from fragment ends
			for(int pos=0; pos<length; ++pos){
				bases[pos].data.distFrom5Prime = pos - softClippedLength[0]; //dist from 5'
				bases[pos].data.distFrom3Prime = fragmentLen - 1 - pos + softClippedLength[0]; //dist from 3'
			}
		}
	} else {
		//treat as single end
		int l = length - 1;
		if(isReverseStrand){
			//not in pair & reverse
			//Hence distance from 3' is just pos
			//And distance from 5' is just len - pos - 1
			for(int pos=0; pos<length; ++pos){
				bases[pos].data.distFrom5Prime = l - pos - softClippedLength[1]; //dist from 5'
				bases[pos].data.distFrom3Prime = pos - softClippedLength[0]; //dist from 3'
			}
		} else {
			//not in pair & forward
			//Hence distance from 5' is just pos
			//And distance from 3' is given by len - pos - 1
			for(int pos=0; pos<length; ++pos){
				bases[pos].data.distFrom5Prime = pos - softClippedLength[0]; //dist from 5'
				bases[pos].data.distFrom3Prime = l - pos - softClippedLength[1]; //dist from 3'
			}
		}
	}
};

void TAlignment::parse(TGenotypeMap & genoMap, TQualityMap & qualityMap){
	if(!parsed){
		if(!storageInitialized)
			throw "Alignment storage was not initialized!";

		//first parse bases and qualities
		parseBasesQualities(genoMap, qualityMap);

		//then update distances from ends
		setDistancesFromEnds();

		//fill context for each base
		fillContext(genoMap);

		//set maping quality and whether read is first or second
		for(int d=0; d<length; ++d){
			bases[d].data.mappingQuality = mappingQuality;
			bases[d].data.fragmentLength = fragmentLength;
			bases[d].setSecondMate(isSecondMate);
			bases[d].setReverseStrand(isReverseStrand);
		}

		parsed = true;
	}
};

void TAlignment::copyDataToBase(TBase & base, const char baseAsChar, const char qualAsChar, TGenotypeMap & genoMap, TQualityMap & qualityMap){
	base.data.base = genoMap.getBase(baseAsChar);
	base.data.qualityAsPhredInt = qualityMap.qualityToPhredInt(qualAsChar);
	base.errorRate = qualityMap.phredIntToError(base.data.qualityAsPhredInt);
};

void TAlignment::parseBasesQualities(TGenotypeMap & genoMap, TQualityMap & qualityMap){
	// iterate over CigarOps
	int d = 0; //index regarding data structures and inside read
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
				for(unsigned int i=0; i<op.Length; ++i, ++d, ++p){
					copyDataToBase(bases[d], bamAlignment.QueryBases[d], bamAlignment.Qualities[d], genoMap, qualityMap);
					bases[d].setAligned(true);
					bases[d].alignedPos = p;
				}
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
				break;

			//for 'S' - soft clip: ignore bases, but increase k
			case (BamTools::Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
				//add bases to softclipped entries
				for(unsigned int i=0; i<op.Length; ++i, ++d){
					//soft-clipped bases on 5' are before bamAlignment.Position
					softClippedBase[softClippedEntry][softClippedLength[softClippedEntry]] = bamAlignment.QueryBases[d];
					softClippedQuality[softClippedEntry][softClippedLength[softClippedEntry]] = bamAlignment.Qualities[d];
					++softClippedLength[softClippedEntry];

					//need to initialize quality for quality filter and bases for context
					copyDataToBase(bases[d], bamAlignment.QueryBases[d], bamAlignment.Qualities[d], genoMap, qualityMap);					bases[d].setAligned(false);
					bases[d].alignedPos = -1;
				}
				break;

			//for 'I' - insertion: copy bases, but put aligned pos to -1
			case (BamTools::Constants::BAM_CIGAR_INS_CHAR)      :
				for(unsigned int i=0; i<op.Length; ++i, ++d){
					copyDataToBase(bases[d], bamAlignment.QueryBases[d], bamAlignment.Qualities[d], genoMap, qualityMap);					bases[d].setAligned(false);
					bases[d].alignedPos = -1;
					++numInsertions;
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
				for(unsigned int i=0; i<op.Length; ++i, ++d, ++p){
					copyDataToBase(bases[d], bamAlignment.QueryBases[d], bamAlignment.Qualities[d], genoMap, qualityMap);					bases[d].setAligned(false);
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
	length = d;
	//TODO: check if we need lastPositionPlusOne and lastAlignedPos
	lastAlignedPositionWithRespectToRef = position + p - 1;
	lastPositionPlusOne = position + length;
	lastAlignedPos = p - 1; //why -1? -> same reason as above
	if(length != bamAlignment.Length)
		throw "The lengths of the alignment and the quality scores of read '" + bamAlignment.Name + "' do not match!";

	//calculate relevant fragment length
	if(isProperPair){
		fragmentLength = abs(insertSize) + numInsertions - numDeletions;
	} else {
		fragmentLength = length - softClippedLength[0] - softClippedLength[1];
	}
};

void TAlignment::fillContext(TGenotypeMap & genoMap){
	if(isReverseStrand){
		//reverse
		for(int d=0; d<(length-1); ++d){
			bases[d].data.context = genoMap.contextMapFlipped[bases[d+1].data.base][bases[d].data.base];
		}
		bases[length-1].data.context = genoMap.contextMapFlipped[N][bases[length-1].data.base];
	} else {
		//forward
		bases[0].data.context = genoMap.contextMap[N][bases[0].data.base];
		for(int d=1; d<length; ++d)
			bases[d].data.context = genoMap.contextMap[bases[d-1].data.base][bases[d].data.base];
	}
};

void TAlignment::fillPmdProbabilities(TPMDDoubleStrand* pmdObjects){
	//probabilities of observing changes in BAM file
	//Note: distances are as in original fragment (not BAM file), i.e. in direction of sequencing
	if(!isReverseStrand){ //is forward
		for(int d=0; d<length; ++d){
			bases[d].PMD_CT = pmdObjects[readGroupId].getProbFivePrime(bases[d].data.distFrom5Prime);
			bases[d].PMD_GA = pmdObjects[readGroupId].getProbThreePrime(bases[d].data.distFrom3Prime);
		}
	} else { //is reverse
		for(int d=0; d<length; ++d){
			bases[d].PMD_CT = pmdObjects[readGroupId].getProbThreePrime(bases[d].data.distFrom3Prime);
			bases[d].PMD_GA = pmdObjects[readGroupId].getProbFivePrime(bases[d].data.distFrom5Prime);
		}
	}
};

//--------------------------------------
//functions to modify data
//--------------------------------------
void TAlignment::filterForBaseQualityAsPhredInt(int & minPhredInt, int & maxPhredInt){
	//set base to N if outside quality filter
	for(int d=0; d<length; ++d){
		if(bases[d].data.qualityAsPhredInt < minPhredInt || bases[d].data.qualityAsPhredInt > maxPhredInt){
			bases[d].data.base = N;
		}
	}
	changed = true;
};

void TAlignment::filterForContext(std::map<BaseContext,int> ignoreTheseContexts){
	//set base to N if outside quality filter
	for(int d=0; d<length; ++d){
		if(ignoreTheseContexts.find(bases[d].data.context) != ignoreTheseContexts.end()){
			bases[d].data.base = N;
			bases[d].errorRate = 1.0;
		}
	}
	changed = true;
};

void TAlignment::filterForPrintingBaseQuality(std::string & bases, std::string & qual, int minQualForPrinting, int maxQualForPrinting){
	//set base to N if outside quality filter
	std::string::iterator stringItBase = bases.begin();
	for(std::string::iterator stringIt = qual.begin() ; stringIt < qual.end(); ++stringIt, ++stringItBase){
		if((int) *stringIt < minQualForPrinting){
			*stringIt = (char) minQualForPrinting;
//			*stringItBase = 'N';
		} else if((int) *stringIt > maxQualForPrinting){
			*stringIt = (char) maxQualForPrinting;
//			*stringItBase = 'N';
		}
	}
};

void TAlignment::trimRead(int & trimmingLength3Prime, int & trimmingLength5Prime){
	//TODO: should we treat paired-end data differently, i.e. should be consider the 3' end of the fragment?
	//set base to N at ends of read
	if(isReverseStrand){
		//distance from 3' is just pos
		//distance from 5' is len - pos - 1
		for(int d=0; d<trimmingLength3Prime; ++d)
			bases[d].data.base = N;
		for(int d=0; d<trimmingLength5Prime; ++d)
			bases[length - 1 - d].data.base = N;
	} else {
		//distance from 3' is len - pos - 1
		//distance from 5' is just pos
		for(int d=0; d<trimmingLength3Prime; ++d)
			bases[length - 1 - d].data.base = N;
		for(int d=0; d<trimmingLength5Prime; ++d)
			bases[d].data.base = N;
	}
	changed = true;
};

void TAlignment::fillReadGroupInfo(int & readGroupID){
	for(int d=0; d<length; ++d){
		bases[d].data.readGroup = readGroupID;
	}
};

void TAlignment::binQualityScores(TQualityMap & qualityMap){
	//make sure read is parsed
	if(!parsed) throw "Read was not parsed!";

	//bin quality scores as done by Illumina
	for(int d=0; d<length; ++d){
		bases[d].errorRate = qualityMap.phredIntToIlluminaError(bases[d].data.qualityAsPhredInt);
	}
	changed = true;
};

void TAlignment::updateOptionalSamField(std::string tag, float value){
	if(bamAlignment.HasTag(tag) == false) bamAlignment.AddTag(tag, "f", value);
	else bamAlignment.EditTag(tag, "f", value);
	changed = true;
};

void TAlignment::updateOptionalSamField(std::string tag, std::string value){
	if(bamAlignment.HasTag(tag) == false) bamAlignment.AddTag(tag, "Z", value);
	else bamAlignment.EditTag(tag, "Z", value);
};

void TAlignment::downsampleAlignment(double& fraction, TRandomGenerator& randomGenerator, TQualityMap & qualMap){
	for(int d=0; d<length; ++d){
		double r = randomGenerator.getRand();
		if(r < fraction){
			bases[d].data.base = N;
			bases[d].errorRate = qualMap.qualityToError(0);
		}
	}
	changed = true;
};

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
	if(!isReverseStrand){ //forward
		for(int d=0; d<length; ++d){
			if(bases[d].isAligned() && bases[d].data.base != N){
				ref = genoMap.getBase(referenceSequence[bases[d].alignedPos]);
				pmdTables.addFromFivePrime(readGroupId, bases[d].data.distFrom5Prime, ref, bases[d].data.base);
				pmdTables.addFromThreePrime(readGroupId, bases[d].data.distFrom3Prime, ref, bases[d].data.base);
			}
		}
	} else { //reverse
		for(int d=0; d<length; ++d){
			if(bases[d].isAligned() && bases[d].data.base != N){
				ref = genoMap.flipBase(referenceSequence[bases[d].alignedPos]);
				read = genoMap.baseToFlippedBase[bases[d].data.base];
				pmdTables.addFromFivePrime(readGroupId, bases[d].data.distFrom5Prime, ref, read);
				pmdTables.addFromThreePrime(readGroupId, bases[d].data.distFrom3Prime, ref, read);
			}
		}
	}};

void TAlignment::recalibrateWithPMD(TRecalibration* recalObject, TQualityMap & qualMap){
	//make sure read is parsed and has reference
	if(!parsed) throw "Read was not parsed!";
	if(!hasReference) throw "Reference was not added!";

	for(int d=0; d<length; ++d){
		if(bases[d].isAligned() && bases[d].data.base != N){
			//alignment is already recalibrated, just need to add effect of PMD
			if(bases[d].data.base == T && referenceSequence[bases[d].alignedPos] == 'C')
				bases[d].errorRate = 1.0 - ((1.0 - bases[d].errorRate)*(1.0 - bases[d].PMD_CT)); //this is mapDamage2, Krishna: qual*(1-pmdCT) + (1-qual)*pmdCT;
			else if(bases[d].data.base == A && referenceSequence[bases[d].alignedPos] == 'G')
				bases[d].errorRate = 1.0 - ((1.0 - bases[d].errorRate)*(1.0 - bases[d].PMD_GA)); //this is mapDamage2, Krishna: qual*(1-pmdGA) + (1-qual)*pmdGA;
		} else {
			bases[d].errorRate = qualMap.phredIntToError(bases[d].data.qualityAsPhredInt);
		}
	}

	recalibrated = true;
	changed = true;
};

double TAlignment::calculatePMDS(double & pi, TPMDDoubleStrand* pmdObjects){
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
		if(bases[d].isAligned() && bases[d].data.base != N && referenceSequence[bases[d].alignedPos] != 'N'){
			//Prepare variables
			epsThird = bases[d].errorRate / 3.0;
			fourEpsThird = 4.0 * epsThird;

			//calc likelihoods
			if(referenceSequence[bases[d].alignedPos] == 'A'){
				if(bases[d].data.base == A){
					probPMD = 1.0 - bases[d].errorRate - pi + fourEpsThird*pi + bases[d].PMD_GA*pi/3.0*(1.0-fourEpsThird);
					probNoPMD = 1.0 - bases[d].errorRate - pi + fourEpsThird*pi;
				}
				else if(bases[d].data.base == C){ //ok
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi - pi*bases[d].PMD_CT*(fourEpsThird-1.0);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else if(bases[d].data.base == G){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + pi*bases[d].PMD_GA*(fourEpsThird-1.0);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else { //T
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + pi*bases[d].PMD_CT*(1.0-fourEpsThird);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
			}
			else if (referenceSequence[bases[d].alignedPos] == 'C'){
				if(bases[d].data.base == A){
					probPMD = bases[d].errorRate + pi - 2.0*bases[d].errorRate*pi + pi*bases[d].PMD_GA*(1.0-fourEpsThird);
					probNoPMD = bases[d].errorRate + pi - 2.0*bases[d].errorRate*pi;
				}
				else if(bases[d].data.base == C){
					probPMD = 1.0 - pi - bases[d].errorRate + fourEpsThird*pi + (1.0-pi)*bases[d].PMD_CT*(fourEpsThird-1.0);
					probNoPMD = 1.0 - pi - bases[d].errorRate + fourEpsThird*pi;
				}
				else if(bases[d].data.base == G){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + bases[d].PMD_GA*(fourEpsThird*pi - pi);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else { //T
					probPMD = epsThird + (1.0-pi)*bases[d].PMD_CT*(1.0-fourEpsThird);
					probNoPMD = epsThird;
				}
			}
			else if (referenceSequence[bases[d].alignedPos] == 'G'){
				if(bases[d].data.base == A){
					probPMD = bases[d].PMD_GA*(3.0-3.0*pi+4.0*bases[d].errorRate+4.0*bases[d].errorRate*pi) + bases[d].errorRate - fourEpsThird*pi + pi;
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else if(bases[d].data.base == C){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + pi*bases[d].PMD_CT*(fourEpsThird - 1.0);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else if(bases[d].data.base == G){
					probPMD = 1.0 - bases[d].errorRate - pi + fourEpsThird*pi + (1.0-pi)*bases[d].PMD_GA*(fourEpsThird-1.0);
					probNoPMD = 1.0 - bases[d].errorRate - pi + fourEpsThird*pi;
				}
				else { //T
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + pi*bases[d].PMD_CT*(1.0-fourEpsThird);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
			}
			else { //T
				if(bases[d].data.base == A){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi - epsThird*pi*bases[d].PMD_CT + pi*bases[d].PMD_GA*(1.0-bases[d].errorRate);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else if(bases[d].data.base == C){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + pi*bases[d].PMD_CT*(fourEpsThird - 1.0);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else if(bases[d].data.base == G){
					probPMD = bases[d].errorRate - fourEpsThird*pi + pi + pi*bases[d].PMD_GA*(fourEpsThird - 1.0);
					probNoPMD = bases[d].errorRate - fourEpsThird*pi + pi;
				}
				else { //T
					probPMD = 1.0 - bases[d].errorRate - pi + fourEpsThird*pi + bases[d].PMD_CT*(pi/3.0-fourEpsThird*pi/3.0);
					probNoPMD = 1.0 - bases[d].errorRate - pi + fourEpsThird*pi;
				}
			}

			//now add to PMDS
			PMDS += log(probPMD/probNoPMD);
		}
	}

	return PMDS;
};

void TAlignment::removeSoftClippedBases(TSoftClippingData & softClippingData){
	//check if there is softclipping
	if(softClippingData.softClippingLength > 0){
		//adapt CIGAR string
		std::vector<BamTools::CigarOp>::const_iterator cigarIter = bamAlignment.CigarData.begin();
		std::vector<BamTools::CigarOp>::const_iterator cigarEnd  = bamAlignment.CigarData.end();

		bamAlignment.CigarData.clear();
		for(; cigarIter != cigarEnd; ++cigarIter ){
			const BamTools::CigarOp& op = (*cigarIter);
			switch ( op.Type ) {

				// for 'M', '=' or 'X': just copy
				case (BamTools::Constants::BAM_CIGAR_MATCH_CHAR)    :
						bamAlignment.CigarData.push_back(BamTools::CigarOp(BamTools::Constants::BAM_CIGAR_MATCH_CHAR, op.Length));
						break;
				case (BamTools::Constants::BAM_CIGAR_SEQMATCH_CHAR) :
						bamAlignment.CigarData.push_back(BamTools::CigarOp(BamTools::Constants::BAM_CIGAR_SEQMATCH_CHAR, op.Length));
						break;
				case (BamTools::Constants::BAM_CIGAR_MISMATCH_CHAR) :
						bamAlignment.CigarData.push_back(BamTools::CigarOp(BamTools::Constants::BAM_CIGAR_MISMATCH_CHAR, op.Length));
						break;
				case (BamTools::Constants::BAM_CIGAR_INS_CHAR)      :
						bamAlignment.CigarData.push_back(BamTools::CigarOp(BamTools::Constants::BAM_CIGAR_INS_CHAR, op.Length));
						break;
				case (BamTools::Constants::BAM_CIGAR_DEL_CHAR) :
						bamAlignment.CigarData.push_back(BamTools::CigarOp(BamTools::Constants::BAM_CIGAR_DEL_CHAR, op.Length));
						break;
				case (BamTools::Constants::BAM_CIGAR_REFSKIP_CHAR) :
						bamAlignment.CigarData.push_back(BamTools::CigarOp(BamTools::Constants::BAM_CIGAR_REFSKIP_CHAR, op.Length));
						break;
				// for 'S': ignore
				case (BamTools::Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
						break;
				// for 'H' - hard clip: do nothing as these bases are not present in SEQ
				case (BamTools::Constants::BAM_CIGAR_HARDCLIP_CHAR) :
						bamAlignment.CigarData.push_back(BamTools::CigarOp(BamTools::Constants::BAM_CIGAR_HARDCLIP_CHAR, op.Length));
						break;;
				// invalid CIGAR op-code
				default:
					throw (std::string) "CIGAR operation type '" + op.Type + "' not supported!";
			}
		}

		bamAlignment.QueryBases = softClippingData.middleBases;
		bamAlignment.Qualities = softClippingData.middleQualities;
	}

	changed = true;
};

int TAlignment::measureOverlap(){
	//make sure read is parsed
	if(!parsed) throw "Read was not parsed!";

	if(isProperPair){
		if(!isReverseStrand){
			int k = length - softClippedLength[1] - 1;
			while(bases[k].alignedPos < 0)
				--k;
			int endPos = position + bases[k].alignedPos;
			int overlap = endPos - matePosition;

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
};

int TAlignment::getUsableLength(const int minQual, const int maxQual){
	int counter = bamAlignment.Length;
	for(int d=0; d<bamAlignment.Length; ++d){
		if(bamAlignment.QueryBases.at(d) == N || bamAlignment.Qualities.at(d) < minQual || bamAlignment.Qualities.at(d) > maxQual){
			counter -= 1;
		}
	}

	return counter;
}

void TAlignment::addToQualityTable(TQualityTable & qualTable, TQualityMap & qualMap){
	//make sure read is parsed
	if(!parsed) throw "Read was not parsed!";
	for(int d=0; d<length; ++d){
		const int qual = qualMap.errorToQuality(bases[d].errorRate);
		qualTable.add(qual);
	}
};

void TAlignment::addToContextStats(TContextStats & contextStats){
	for(int d=0; d<length; ++d){
		contextStats.add(bases[d]);
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

void TAlignment::setIsProperPair(const bool & ok){
	bamAlignment.SetIsProperPair(ok);
	isProperPair = false;
}

void TAlignment::save(BamTools::BamWriter & bamWriter, TGenotypeMap & genoMap, int minQualForPrinting, int maxQualForPrinting, TQualityMap & qualMap){
	if(changed){
		//means that read has been modified.
		//Currently quality recalibration and PMDS flag is the only possible change.
		//But will need to think how to deal with merging and such...

		//assume that only bases and quality scores where changed
		std::string tmpString, tmpString2;
		for(int d=0; d<length; ++d){
			tmpString += genoMap.baseToChar[bases[d].data.base];
			tmpString2 += (char) qualMap.errorToQuality(bases[d].errorRate);
		}
		bamAlignment.QueryBases = tmpString;
		bamAlignment.Qualities = tmpString2;
	}

	//make sure quality are printed within range
	filterForPrintingBaseQuality(bamAlignment.QueryBases, bamAlignment.Qualities, minQualForPrinting, maxQualForPrinting);

	//now write alignment
	if(!bamWriter.SaveAlignment(bamAlignment))
		throw "Read '" + bamAlignment.Name + "' could not be written!";
};

void TAlignment::print(TGenotypeMap & genoMap, TQualityMap & qualMap){
	std::cout << "NAME:\t" << alignmentName << std::endl;
	std::cout << "aligned bases " << bamAlignment.AlignedBases << std::endl;
	std::cout << "LEN:\t" << length << std::endl;

	//print bases
	std::cout << "SEQ:\t";
	for(int d=0; d<length; ++d)
		std::cout << genoMap.getBaseAsChar(bases[d].data.base);
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
		if(bases[d].isAligned())
			std::cout << bases[d].alignedPos;
		else
			std::cout << "-";
	}
	std::cout << std::endl;

	//print dist from 3'
	std::cout << "dist 3':\t";
	for(int d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		std::cout << bases[d].data.distFrom3Prime;
	}
	std::cout << std::endl;

	//print dist from 5'
	std::cout << "dist 5':\t";
	for(int d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		std::cout << bases[d].data.distFrom5Prime;
	}
	std::cout << std::endl;
}
