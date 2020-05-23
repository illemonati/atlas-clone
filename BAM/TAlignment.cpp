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
	refID = 0;
	readGroupID = 0;
	position = 0;
	matePosition = 0;
	lastAlignedPos = 0;
	lastAlignedPositionWithRespectToRef = 0;
	insertSize_TLEN = 0;
	isReverseStrand = false;
	isPaired = false;
	isProperPair = false;
	mappingQuality = 0;
	isSecondMate = false;
	parsed = false;
	changed = false;
	storageInitialized = false;
	std::string referenceSequence = "";
	hasReference = false;

	bases = nullptr;
	alignedPosition = nullptr;
	numInsertions = -1;
	numDeletions = -1;
	softClippedEntry = 0;
	softClippedLength = NULL;
	softClippedBase = NULL;
	softClippedQuality = NULL;
};

TAlignment::TAlignment(uint16_t MaxSize){
	TAlignment();
	maxSize = MaxSize;
	_initStorage();
};

TAlignment::TAlignment(const TAlignment & Alignment){
	TAlignment();

	//initialize storage
	maxSize = Alignment.maxSize;
	_initStorage();
	storageInitialized = true;

	//copy content
	name = Alignment.name;
	length = Alignment.length;
	refID = Alignment.refID;
	readGroupID = Alignment.readGroupID;
	position = Alignment.position;
	lastAlignedPos = Alignment.lastAlignedPos;
	lastAlignedPositionWithRespectToRef = Alignment.lastAlignedPositionWithRespectToRef;
	mappingQuality = Alignment.mappingQuality;
	parsed = Alignment.parsed;
	changed = Alignment.changed;
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
	name.clear();
	length = 0;
	empty = true;
	parsed = false;
	changed = false;
};

void TAlignment::_initStorage(){
	clear();

	//data
	bases = new TBase[maxSize];
	alignedPosition = new uint16_t[maxSize];

	//soft clipped data
	softClippedLength = new int[2];
	softClippedBase = new char*[2];
	softClippedQuality = new char*[2];
	softClippedBase[0] = new char[maxSize];
	softClippedBase[1] = new char[maxSize];
	softClippedQuality[0] = new char[maxSize];
	softClippedQuality[1] = new char[maxSize];

	storageInitialized = true;
};

void TAlignment::_freeStorage(){
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


//--------------------------------------
// functions to fill alignment
//--------------------------------------

//function used by TBamFile to fill alignment
void TAlignment::fill(const	std::string Name,
		  const BAM::TSamFlags Flags,
		  const uint32_t RefID,
		  const uint32_t Position,
		  const uint16_t MappingQuality,
		  const BAM::TCigar Cigar,
		  const uint32_t MateRefID,
		  const int32_t MatePosition,
		  const int32_t InsertSize_TLEN,
		  const std::string Sequence,
		  const std::string Qualities,
		  const uint16_t ReadGroupId){

	//empty alignment
	clear();

	//copy data
	name = Name;
	flags = Flags;
	refID = RefID;
	position = Position;
	mappingQuality = MappingQuality;
	cigar = Cigar;
	mateRefID = MateRefID;
	matePosition = MatePosition;
	insertSize_TLEN = InsertSize_TLEN;
	sequence = Sequence;
	qualities = Qualities;
	readGroupID = ReadGroupId;
	empty = false;
};


void TAlignment::parse(const TGenotypeMap & genoMap, const TQualityMap & qualityMap, const GenotypeLikelihoods::TSequencingErrorModels & seqErrorModels){
	if(!parsed){
		if(!storageInitialized)
			throw "Alignment storage was not initialized!";

		//first parse bases and qualities
		_parseBasesQualities(genoMap, qualityMap);

		//then update distances from ends
		_setDistancesFromEnds();

		//fill context for each base
		_fillContext();

		//set mapping quality and whether read is first or second
		for(int d=0; d<length; ++d){
			bases[d].readGroupID = readGroupID;
			bases[d].mappingQuality = mappingQuality;
			bases[d].fragmentLength = fragmentLength;
			bases[d].setSecondMate(flags.isSecondMate());
			bases[d].setReverseStrand(flags.isReverseStrand());
		}

		//recalibrate
		seqErrorModels.recalibrate(bases, length);

		parsed = true;
		changed = seqErrorModels.recalibrationChangesQualities();
	}
};

void TAlignment::_parseBasesQualities(const TGenotypeMap & genoMap, const TQualityMap & qualityMap){
	//initialize
	int d = 0; //index regarding data structures and inside read
	int p = 0; //index regarding reference position (!= k for indels)
	softClippedEntry = 0; //softclipped bases to be added
	softClippedLength[0] = 0; softClippedLength[1] = 0;
	numInsertions = 0;
	numDeletions = 0;

	//loop over cigar operations
	for(auto& cigarIter : cigar){
		switch ( cigarIter.type ) {

			// for 'M', '=' or 'X': just copy
			case ('M') :
			case ('=') :
			case ('X') :

				//soft-clipped bases on 5' are before bamAlignment.Position
				for(unsigned int i=0; i<cigarIter.length; ++i, ++d, ++p){
					bases[d].base = genoMap.getBase(sequence[d]);
					bases[d].originalQuality_phredInt = qualityMap.qualityToPhredInt(qualities[d]);
					bases[d].setAligned(true);
					alignedPosition[d] = p;
				}
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
				break;

			//for 'S' - soft clip: ignore bases, but increase k
			case ('S') :
				//add bases to softclipped entries
				for(unsigned int i=0; i<cigarIter.length; ++i, ++d){
					//soft-clipped bases on 5' are before bamAlignment.Position
					//TODO: find better way than copying bases with indexes. use string?
					softClippedBase[softClippedEntry][softClippedLength[softClippedEntry]] = sequence[d];
					softClippedQuality[softClippedEntry][softClippedLength[softClippedEntry]] = qualities[d];
					++softClippedLength[softClippedEntry];

					//need to initialize quality for quality filter and bases for context
					bases[d].base = genoMap.getBase(sequence[d]);
					bases[d].originalQuality_phredInt = qualityMap.qualityToPhredInt(qualities[d]);
					bases[d].setAligned(false);
					alignedPosition[d] = -1;
				}
				break;

			//for 'I' - insertion: copy bases, but put aligned pos to -1
			case ('I')      :
				for(unsigned int i=0; i<cigarIter.length; ++i, ++d){
					bases[d].base = genoMap.getBase(sequence[d]);
					bases[d].originalQuality_phredInt = qualityMap.qualityToPhredInt(qualities[d]);
					bases[d].setAligned(false);
					alignedPosition[d] = -1;
					++numInsertions;
				}
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
				break;


			// for 'D' - deletion: just add to position
			case ('D') :
				p += cigarIter.length;
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
				++numDeletions;
				break;

			// for 'N' - skipped region in reference: only advance reference position
			case ('N') :
				p += cigarIter.length;
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
				break;

			// for 'H' or 'P' - hard clip: do nothing as these bases are not present in SEQ
			case ('H') :
			case ('P') :
				break;

			// invalid CIGAR op-code
			default:
				throw (std::string) "CIGAR operation '" + cigarIter.type + "' not supported!";
		}
	}

	//update length and last aligned position
	length = d;
	lastAlignedPositionWithRespectToRef = position + p - 1;
	lastAlignedPos = p - 1; //why -1? -> same reason as above
	if(length != sequence.length())
		throw "The lengths of the alignment and the quality scores of read '" + name + "' do not match!";

	//calculate relevant fragment length
	if(flags.isProperPair()){
		fragmentLength = abs(insertSize_TLEN) + numInsertions - numDeletions;
	} else {
		fragmentLength = length - softClippedLength[0] - softClippedLength[1];
	}
};

void TAlignment::_setDistancesFromEnds(){
	//Set distances in ORIGINAL FRAGMENT (i.e. 5' end is where sequencing started, NOT how it aligns to reference)
	//is it paired-end?
	if(flags.isProperPair()){
		int fragmentLen = abs(insertSize_TLEN) + numInsertions - numDeletions;
		if(flags.isReverseStrand()){
			//reverse (can be either first or second mate, but it's the one that comes second in bam file)
			//and distance from 5' is given as f(end of fragment) = f(len - pos - 1)
			//hence distance from 3' is given by f(dist since beginning of fragment) = f(insert - len + pos)
			int k = abs(fragmentLen) - (length - softClippedLength[1]);
			int l = length - 1 - softClippedLength[1];
			for(int pos=0; pos<length; ++pos){
				bases[pos].distFrom5Prime = l - pos; //dist from 5'
				bases[pos].distFrom3Prime = k + pos; //dist from 3'
			}
		} else {
			//forward (can be either first or second mate, but it's the one that comes first in bam file)
			//Hence distance from 5' is given as a function of pos
			//And distance from 3' is given by (length of fragment) - pos -1
			//NOTE! we ignore indels when calculating distance from 5' since we can not know this info.
			//Luckily, this has only minimal effect since these distances are far from fragment ends
			for(int pos=0; pos<length; ++pos){
				bases[pos].distFrom5Prime = pos - softClippedLength[0]; //dist from 5'
				bases[pos].distFrom3Prime = fragmentLen - 1 - pos + softClippedLength[0]; //dist from 3'
			}
		}
	} else {
		//treat as single end
		int l = length - 1;
		if(flags.isReverseStrand()){
			//not in pair & reverse
			//Hence distance from 3' is just pos
			//And distance from 5' is just len - pos - 1
			for(int pos=0; pos<length; ++pos){
				bases[pos].distFrom5Prime = l - pos - softClippedLength[1]; //dist from 5'
				bases[pos].distFrom3Prime = pos - softClippedLength[0]; //dist from 3'
			}
		} else {
			//not in pair & forward
			//Hence distance from 5' is just pos
			//And distance from 3' is given by len - pos - 1
			for(int pos=0; pos<length; ++pos){
				bases[pos].distFrom5Prime = pos - softClippedLength[0]; //dist from 5'
				bases[pos].distFrom3Prime = l - pos - softClippedLength[1]; //dist from 3'
			}
		}
	}
};

void TAlignment::_fillContext(){
	if(flags.isReverseStrand()){
		//reverse
		for(int d=0; d<(length-1); ++d){
			bases[d].context = bases[d+1].base;
		}
		bases[length-1].context = N;
	} else {
		//forward
		bases[0].context = N;
		for(int d=1; d<length; ++d)
			bases[d].context = bases[d-1].base;
	}
};

//--------------------------------------
//getters
//--------------------------------------
void TAlignment::_updateSequenceAndQualities(const TGenotypeMap & genoMap, const TQualityMap & qualMap){
	if(sequenceAndQualitiesChanged){
		//update according to what is stored in bases
		sequence.clear();
		qualities.clear();

		for(int d=0; d<length; ++d){
			sequence += genoMap.baseToChar[bases[d].base];
			qualities += (char) qualMap.phredIntToQuality(bases[d].recalibratedQualityAsPhredInt);
		}

		sequenceAndQualitiesChanged = false;
	}
}

std::string TAlignment::getSequence(const TGenotypeMap & genoMap, const TQualityMap & qualMap){
	_updateSequenceAndQualities(genoMap, qualMap);
	return sequence;
};

std::string TAlignment::getQualities(const TGenotypeMap & genoMap, const TQualityMap & qualMap){
	_updateSequenceAndQualities(genoMap, qualMap);
	return qualities;
}

//--------------------------------------
//functions to modify data
//--------------------------------------
void TAlignment::filterForBaseQualityAsPhredInt(int & minPhredInt, int & maxPhredInt){
	//set base to N if outside quality filter
	for(int d=0; d<length; ++d){
		if(bases[d].originalQuality_phredInt < minPhredInt || bases[d].originalQuality_phredInt > maxPhredInt){
			bases[d].base = N;
		}
	}
	changed = true;
};

void TAlignment::filterForContext(std::map<BaseContext,int> ignoreTheseContexts){
	//set base to N if outside quality filter

	//TODO: discuss what to do!
	/*
	for(int d=0; d<length; ++d){
		if(ignoreTheseContexts.find(bases[d].context) != ignoreTheseContexts.end()){
			bases[d].base = N;
			bases[d].recalibratedQualityAsPhredInt = 0;
		}
	}
	*/
	changed = true;
};

void TAlignment::_trimRead(int & trimmingLength3Prime, int & trimmingLength5Prime){
	//TODO: should we treat paired-end data differently, i.e. should be consider the 3' end of the fragment?
	//set base to N at ends of read
	if(flags.isReverseStrand()){
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

void TAlignment::addReference(TFastaBuffer & fasta){
	fasta.fill(refID, position, lastAlignedPositionWithRespectToRef, referenceSequence);
	hasReference = true;
};

void TAlignment::binQualityScores(TQualityMap & qualityMap){
	//make sure read is parsed
	if(!parsed) throw "Read was not parsed!";

	//bin quality scores as done by Illumina
	for(int d=0; d<length; ++d){
		bases[d].recalibratedQualityAsPhredInt = qualityMap.phredIntToIlluminaError(bases[d].originalQuality_phredInt);
	}
	changed = true;
};

/*
void TAlignment::updateOptionalSamField(std::string tag, float value){
	if(bamAlignment.HasTag(tag) == false) bamAlignment.AddTag(tag, "f", value);
	else bamAlignment.EditTag(tag, "f", value);
	changed = true;
};

void TAlignment::updateOptionalSamField(std::string tag, std::string value){
	if(bamAlignment.HasTag(tag) == false) bamAlignment.AddTag(tag, "Z", value);
	else bamAlignment.EditTag(tag, "Z", value);
};
*/

void TAlignment::downsampleAlignment(double& fraction, TRandomGenerator& randomGenerator, TQualityMap & qualMap){
	for(int d=0; d<length; ++d){
		double r = randomGenerator.getRand();
		if(r < fraction){
			bases[d].base = N;
			bases[d].recalibratedQualityAsPhredInt = 0;
		}
	}
	changed = true;
};

//--------------------------------------
//functions to fill other classes
//--------------------------------------
void TAlignment::addToPMDTables(GenotypeLikelihoods::TPMDTables & pmdTables, TGenotypeMap & genoMap){
	//make sure read is parsed and has reference
	if(!parsed) throw "Read was not parsed!";
	if(!hasReference) throw "Reference was not added!";

	//tmp variables
	Base ref, read;

	//check if it is forward or reverse strand!
	if(!isReverseStrand){ //forward
		for(int d=0; d<length; ++d){
			if(bases[d].isAligned() && bases[d].base != N){
				ref = genoMap.getBase(referenceSequence[alignedPosition[d]]);
				pmdTables.addFromFivePrime(readGroupID, bases[d].distFrom5Prime, ref, bases[d].base);
				pmdTables.addFromThreePrime(readGroupID, bases[d].distFrom3Prime, ref, bases[d].base);
			}
		}
	} else { //reverse
		for(int d=0; d<length; ++d){
			if(bases[d].isAligned() && bases[d].base != N){
				ref = genoMap.flipBase(referenceSequence[alignedPosition[d]]);
				read = genoMap.baseToFlippedBase[bases[d].base];
				pmdTables.addFromFivePrime(readGroupID, bases[d].distFrom5Prime, ref, read);
				pmdTables.addFromThreePrime(readGroupID, bases[d].distFrom3Prime, ref, read);
			}
		}
	}
};

void TAlignment::addSitesToQualityTransformTable(TQualityTransformTables & QTtables){
	for(int i=0; i<length; ++i){
		if(bases[i].base != N){
			QTtables.add(readGroupID, bases[i].originalQuality_phredInt, bases[i].recalibratedQualityAsPhredInt);
		}
	}
};

void TAlignment::addSitesToQualityTransformTable(GenotypeLikelihoods::TSequencingErrorModels & otherSeqErrors, TQualityTransformTables & QTtables){
	for(int i=0; i<length; ++i){
		if(bases[i].base != N){
			QTtables.add(readGroupID, bases[i].recalibratedQualityAsPhredInt, otherSeqErrors.getPhredInt(bases[i]));
		}
	}
};

void TAlignment::recalibrateWithPMD(GenotypeLikelihoods::TGenotypeLikelihoodCalculator & GLCalculator){
	GLCalculator.recalibrateWithPMD(bases, length);
	changed = true;
};

double TAlignment::calculatePMDS(double & pi, GenotypeLikelihoods::TGenotypeLikelihoodCalculator & GLCalculator){
	//make sure read is parsed and has reference
	if(!parsed) throw "Read was not parsed!";
	if(!hasReference) throw "Reference was not added!";

	//calcukate PMDS (is in log)
	double PMDS = 0.0;
	for(int d=0; d<length; ++d){
		//limit to aligned positions
		if(bases[d].isAligned() && bases[d].base != N && referenceSequence[alignedPosition[d]] != 'N'){
			PMDS +=  GLCalculator.calculatePMDS(bases[d], referenceSequence[alignedPosition[d]], pi);
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
			while(!bases[k].isAligned())
				--k;
			int endPos = position + alignedPosition[k];
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
		qualTable.add(bases[d].recalibratedQualityAsPhredInt);
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
void TAlignment::setIsProperPair(const bool & ok){
	bamAlignment.SetIsProperPair(ok);
	isProperPair = false;
}

void TAlignment::save(BamTools::BamWriter & bamWriter, const TGenotypeMap & genoMap, const TQualityMap & qualMap){
	if(changed){
		//means that read has been modified.

		//assume that only bases and quality scores where changed
		std::string tmpString, tmpString2;
		for(int d=0; d<length; ++d){
			tmpString += genoMap.baseToChar[bases[d].base];
			tmpString2 += (char) qualMap.phredIntToQuality(bases[d].recalibratedQualityAsPhredInt);
		}
		bamAlignment.QueryBases = tmpString;
		bamAlignment.Qualities = tmpString2;
	}

	//now write alignment
	if(!bamWriter.SaveAlignment(bamAlignment))
		throw "Read '" + bamAlignment.Name + "' could not be written!";
};

void TAlignment::print(TGenotypeMap & genoMap, TQualityMap & qualMap){
	std::cout << "NAME:\t" << name << std::endl;
	std::cout << "LEN:\t" << length << std::endl;

	//print bases
	std::cout << "SEQ:\t";
	for(int d=0; d<length; ++d)
		std::cout << genoMap.getBaseAsChar(bases[d].base);
	std::cout << std::endl;

	//print qualities
	std::cout << "QUAL:\t";
	for(int d=0; d<length; ++d)
		std::cout << (char) qualMap.phredIntToQuality(bases[d].recalibratedQualityAsPhredInt);
	std::cout << std::endl;

	//print aligned pos
	std::cout << "POS:\t";
	for(int d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		if(bases[d].isAligned())
			std::cout << alignedPosition[d];
		else
			std::cout << "-";
	}
	std::cout << std::endl;

	//print dist from 3'
	std::cout << "dist 3':\t";
	for(int d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		std::cout << bases[d].distFrom3Prime;
	}
	std::cout << std::endl;

	//print dist from 5'
	std::cout << "dist 5':\t";
	for(int d=0; d<length; ++d){
		if(d>0) std::cout << ",";
		std::cout << bases[d].distFrom5Prime;
	}
	std::cout << std::endl;
}
