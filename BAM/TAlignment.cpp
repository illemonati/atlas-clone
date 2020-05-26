/*
 * TAlignment.cpp
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */

#include "TAlignment.h"

namespace BAM{

TAlignment::TAlignment(){
	//details
	empty = true;
	refID = 0;
	readGroupID = 0;
	position = 0;
	matePosition = 0;
	insertSize_TLEN = 0;
	mappingQuality = 0;
	parsed = false;
	changed = false;
	std::string referenceSequence = "";
	hasReference = false;

	lastAlignedPos = 0;
	lastAlignedPositionWithRespectToRef = 0;

	alignedPosition = nullptr;

	/*
	softClippedEntry = 0;
	softClippedLength = NULL;
	softClippedBase = NULL;
	softClippedQuality = NULL;
	*/
};

void TAlignment::clear(){
	position = -1;
	name.clear();
	cigar.clear();
	empty = true;
	parsed = false;
	changed = false;
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

	//set fragment length
	if(flags.isProperPair()){
		fragmentLength = abs(insertSize_TLEN) + cigar.lengthInserted() - cigar.lengthDeleted();
	} else {
		fragmentLength = cigar.lengthSequenced();
	}
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
		for(auto& b : bases){
			b.readGroupID = readGroupID;
			b.mappingQuality = mappingQuality;
			b.fragmentLength = fragmentLength;
			b.setSecondMate(flags.isSecondMate());
			b.setReverseStrand(flags.isReverseStrand());
		}

		//recalibrate
		seqErrorModels.recalibrate(bases);

		parsed = true;
		changed = seqErrorModels.recalibrationChangesQualities();
	}
};

void TAlignment::_parseBasesQualities(const TGenotypeMap & genoMap, const TQualityMap & qualityMap){
	//initialize
	bases.resize(cigar.lengthSequenced());
	alignedPosition.resize(cigar.lengthSequenced());
	int d = 0; //index regarding data structures and inside read
	int p = 0; //index regarding reference position (!= k for indels)
	uint8_t softClippedEntry = 0; //softclipped bases to be added
	softClippedLength[0] = 0; softClippedLength[1] = 0;

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
				softClippedLength[softClippedEntry] += cigarIter.length;
				for(unsigned int i=0; i<cigarIter.length; ++i, ++d){
					//soft-clipped bases on 5' are before bamAlignment.Position
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
				}
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
				break;


			// for 'D' - deletion: just add to position
			case ('D') :
				p += cigarIter.length;
				softClippedEntry = 1; //soft clipped bases can now only occur at the end!
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

	//calculate relevant fragment length

	//update length and last aligned position
	lastAlignedPositionWithRespectToRef = position + p - 1;
	lastAlignedPos = p - 1; //why -1? -> same reason as above
};

void TAlignment::_setDistancesFromEnds(){
	//Set distances in ORIGINAL FRAGMENT (i.e. 5' end is where sequencing started, NOT how it aligns to reference)
	int length = cigar.lengthSequenced();

	//is it paired-end?
	if(flags.isProperPair()){
		if(flags.isReverseStrand()){
			//reverse (can be either first or second mate, but it's the one that comes second in bam file)
			//and distance from 5' is given as f(end of fragment) = f(len - pos - 1)
			//hence distance from 3' is given by f(dist since beginning of fragment) = f(insert - len + pos)
			int k = abs(fragmentLength) - (cigar.lengthSequenced() - softClippedLength[1]);
			int l = cigar.lengthSequenced() - 1 - softClippedLength[1];
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
				bases[pos].distFrom3Prime = fragmentLength - 1 - pos + softClippedLength[0]; //dist from 3'
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
		for(int d=0; d<(cigar.lengthSequenced()-1); ++d){
			bases[d].context = bases[d+1].base;
		}
		bases[cigar.lengthSequenced()-1].context = N;
	} else {
		//forward
		bases[0].context = N;
		for(int d=1; d<cigar.lengthSequenced(); ++d)
			bases[d].context = bases[d-1].base;
	}
};

void TAlignment::addReference(TFastaBuffer & fasta){
	fasta.fill(refID, position, lastAlignedPositionWithRespectToRef, referenceSequence);
	hasReference = true;
};

//--------------------------------------
//getters
//--------------------------------------
void TAlignment::_updateSequenceAndQualities(const TGenotypeMap & genoMap, const TQualityMap & qualMap){
	if(sequenceAndQualitiesChanged || qualMap.qualityLimitSet){
		//update according to what is stored in bases
		sequence.clear();
		qualities.clear();

		for(auto& b : bases){
			sequence += genoMap.baseToChar[b.base];
			qualities += (char) qualMap.phredIntToQuality(b.recalibratedQualityAsPhredInt);
		}

		sequenceAndQualitiesChanged = false;
	}
};

std::string TAlignment::getSequence(const TGenotypeMap & genoMap, const TQualityMap & qualMap){
	_updateSequenceAndQualities(genoMap, qualMap);
	return sequence;
};

std::string TAlignment::getQualities(const TGenotypeMap & genoMap, const TQualityMap & qualMap){
	_updateSequenceAndQualities(genoMap, qualMap);
	return qualities;
};

double TAlignment::calculatePMDS(const double pi, GenotypeLikelihoods::TGenotypeLikelihoodCalculator & GLCalculator, const TGenotypeMap & genoMap){
	//make sure read is parsed and has reference
	if(!parsed) throw "Read was not parsed!";
	if(!hasReference) throw "Reference was not added!";

	//calcukate PMDS (is in log)
	double PMDS = 0.0;
	for(size_t d=0; d<bases.size(); ++d){
		//limit to aligned positions
		if(bases[d].isAligned() && bases[d].base != N && referenceSequence[alignedPosition[d]] != 'N'){
			Base ref = genoMap.getBase(referenceSequence[alignedPosition[d]]);
			PMDS +=  GLCalculator.calculatePMDS(bases[d], ref, pi);
		}
	}

	return PMDS;
};

int TAlignment::measureOverlap(){
	//TODO: check that it works! What is the purpos eof only doing this on reverse reads? Would probably be better to have this work for all reads and then only call it on correct reads?
	//make sure read is parsed
	if(!parsed) throw "Read was not parsed!";

	if(flags.isProperPair()){
		if(!flags.isReverseStrand()){
			int k = cigar.lengthSequenced() - 1;
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

//--------------------------------------------
//filters and other functions to modify data
//--------------------------------------------
void TAlignment::filterForBaseQualityAsPhredInt(const int & minPhredInt, const int & maxPhredInt){
	//set base to N if outside quality filter
	for(auto& b : bases){
		if(b.recalibratedQualityAsPhredInt < minPhredInt || b.recalibratedQualityAsPhredInt > maxPhredInt){
			b.base = N;
			b.recalibratedQualityAsPhredInt = 0;
		}
	}
	changed = true;
};

void TAlignment::filterForContext(const std::map<BaseContext,int> & ignoreTheseContexts, const TGenotypeMap & genoMap){
	//set base to N if outside quality filter
	for(auto& b : bases){
		BaseContext c = genoMap.getContext(b.base, b.context);
		if(ignoreTheseContexts.find(c) != ignoreTheseContexts.end()){
			b.base = N;
			b.recalibratedQualityAsPhredInt = 0;
		}
	}

	changed = true;
};

void TAlignment::trimRead(const int & trimmingLength3Prime, const int & trimmingLength5Prime){
	for(auto& b : bases){
		if(b.distFrom3Prime < trimmingLength3Prime || b.distFrom5Prime < trimmingLength5Prime){
			b.base = N;
			b.recalibratedQualityAsPhredInt = 0;
		}
	}

	changed = true;
};

void TAlignment::removeSoftClippedBases(TSoftClippingData & softClippingData){
	//check if there is softclipping
	if(softClippingData.softClippingLength > 0){
		auto bIter = bases.begin();
		for(auto& cigarIter : cigar){
			if(cigarIter.type == 'S'){
				//remove bases
				bIter = bases.erase(bIter, bIter + cigarIter.length);
			} else {
				//just advance position
				bIter += cigarIter.length;
			}
		}
	}

	changed = true;
};

void TAlignment::binQualityScores(TQualityMap & qualityMap){
	//make sure read is parsed
	if(!parsed) throw "Read was not parsed!";

	//bin quality scores as done by Illumina
	for(auto& b : bases){
		b.recalibratedQualityAsPhredInt = qualityMap.phredIntToIlluminaError(b.recalibratedQualityAsPhredInt);
	}
	changed = true;
};

void TAlignment::recalibrateWithPMD(GenotypeLikelihoods::TGenotypeLikelihoodCalculator & GLCalculator){
	GLCalculator.recalibrateWithPMD(bases);
	changed = true;
};

void TAlignment::setIsProperPair(const bool & ok){
	flags.setIsProperPair(ok);
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

void TAlignment::downsampleAlignment(double& fractionToKeep, TRandomGenerator& randomGenerator, TQualityMap & qualMap){
	for(auto& b : bases){
		double r = randomGenerator.getRand();
		if(r > fractionToKeep){
			b.base = N;
			b.recalibratedQualityAsPhredInt = 0;
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
	if(!flags.isReverseStrand()){ //forward
		int d = 0;
		for(auto& b : bases){
			if(b.isAligned() && b.base != N){
				ref = genoMap.getBase(referenceSequence[alignedPosition[d]]);
				pmdTables.addFromFivePrime(readGroupID, b.distFrom5Prime, ref, b.base);
				pmdTables.addFromThreePrime(readGroupID, b.distFrom3Prime, ref, b.base);
			}
			++d;
		}
	} else { //reverse
		int d = 0;
		for(auto& b : bases){
			if(b.isAligned() && b.base != N){
				ref = genoMap.flipBase(referenceSequence[alignedPosition[d]]);
				read = genoMap.baseToFlippedBase[b.base];
				pmdTables.addFromFivePrime(readGroupID, b.distFrom5Prime, ref, read);
				pmdTables.addFromThreePrime(readGroupID, b.distFrom3Prime, ref, read);
			}
			++d;
		}
	}
};

void TAlignment::addSitesToQualityTransformTable(TQualityTransformTables & QTtables){
	for(auto& b : bases){
		if(b.base != N){
			QTtables.add(readGroupID, b.originalQuality_phredInt, b.recalibratedQualityAsPhredInt);
		}
	}
};

void TAlignment::addSitesToQualityTransformTable(GenotypeLikelihoods::TSequencingErrorModels & otherSeqErrors, TQualityTransformTables & QTtables){
	for(auto& b : bases){
		if(b.base != N){
			QTtables.add(readGroupID, b.recalibratedQualityAsPhredInt, otherSeqErrors.getPhredInt(b));
		}
	}
};

void TAlignment::addToQualityTable(TQualityTable & qualTable, TQualityMap & qualMap){
	//make sure read is parsed
	if(!parsed) throw "Read was not parsed!";
	for(auto& b : bases){
		qualTable.add(b.recalibratedQualityAsPhredInt);
	}
};

void TAlignment::addToContextStats(TContextStats & contextStats){
	for(auto& b : bases){
		contextStats.add(b);
	}
};

//--------------------------------------------
//functions to write / print alignment
//--------------------------------------------
void TAlignment::print(TGenotypeMap & genoMap, TQualityMap & qualMap){
	std::cout << "NAME:\t" << name << std::endl;
	std::cout << "LEN:\t" << bases.size() << std::endl;

	//print bases
	std::cout << "SEQ:\t";
	for(auto& b : bases){
		std::cout << genoMap.getBaseAsChar(b.base);
	}
	std::cout << std::endl;

	//print qualities
	std::cout << "QUAL:\t";
	for(auto& b : bases){
		std::cout << (char) qualMap.phredIntToQuality(b.recalibratedQualityAsPhredInt);
	std::cout << std::endl;

	//print aligned pos
	std::cout << "POS:\t";
	for(size_t d=0; d<bases.size(); ++d){
		if(d>0) std::cout << ",";
		if(bases[d].isAligned())
			std::cout << alignedPosition[d];
		else
			std::cout << "-";
	}
	std::cout << std::endl;

	//print dist from 3'
	std::cout << "dist 3':\t";
	bool first = true;
	for(auto& b : bases){
		if(first){
			first = false;
		} else {
			std::cout << ",";
		}
		std::cout << b.distFrom3Prime;
	}
	std::cout << std::endl;

	//print dist from 5'
	std::cout << "dist 5':\t";
	bool first = true;
	for(auto& b : bases){
		if(first){
			first = false;
		} else {
			std::cout << ",";
		}
		std::cout << b.distFrom5Prime;
	}
	std::cout << std::endl;
};

}; //end namespace
