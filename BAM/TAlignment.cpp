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
	_empty = true;
	_readGroupID = 0;
	_fragmentLength = 0;
	_insertSize_TLEN = 0;
	_mappingQuality = 0;
	_parsed = false;
	_changed = false;
	_sequenceAndQualitiesChanged = false;
	_hasReference = false;
	_lastAlignedPos = 0;
	_lastAlignedPositionWithRespectToRef = 0;
};

void TAlignment::clear(){
	_name.clear();
	_cigar.clear();
	_mateGenomicPosition.clear();
	_lastAlignedPositionWithRespectToRef.clear();
	_empty = true;
	_parsed = false;
	_changed = false;
};

//--------------------------------------
// functions to fill alignment
//--------------------------------------
//function used by TBamFile to fill alignment
void TAlignment::fill(const	std::string Name,
		  const TSamFlags Flags,
		  const uint32_t RefID,
		  const uint32_t Position,
		  const uint16_t MappingQuality,
		  const TCigar Cigar,
		  const uint32_t MateRefID,
		  const uint32_t MatePosition,
		  const int32_t InsertSize_TLEN,
		  const std::string Sequence,
		  const std::string Qualities,
		  const uint16_t ReadGroupId){

	//empty alignment
	clear();

	//copy data
	_name = Name;
	_flags = Flags;
	update(RefID, Position);
	_mappingQuality = MappingQuality;
	_cigar = Cigar;
	_mateGenomicPosition.update(MateRefID, MatePosition);
	_insertSize_TLEN = InsertSize_TLEN;
	_sequence = Sequence;
	_qualities = Qualities;
	_readGroupID = ReadGroupId;
	_empty = false;

	//set fragment length
	if(_flags.isProperPair()){
		_fragmentLength = abs(_insertSize_TLEN) + _cigar.lengthInserted() - _cigar.lengthDeleted();
	} else {
		_fragmentLength = _cigar.lengthSequenced();
	}
};

void TAlignment::parse(const TGenotypeMap & genoMap, const TQualityMap & qualityMap){
	//first parse bases and qualities
	_parseBasesQualities(genoMap, qualityMap);

	//then update distances from ends
	_setDistancesFromEnds();

	//fill context for each base
	_fillContext();

	//set mapping quality and whether read is first or second
	for(auto& b : _bases){
		b.readGroupID = _readGroupID;
		b.mappingQuality = _mappingQuality;
		b.fragmentLength = _fragmentLength;
		b.setSecondMate(_flags.isSecondMate());
		b.setReverseStrand(_flags.isReverseStrand());
	}

	_parsed = true;
	_sequenceAndQualitiesChanged = false;
};


void TAlignment::parse(const TGenotypeMap & genoMap, const TQualityMap & qualityMap, const GenotypeLikelihoods::TSequencingErrorModels & seqErrorModels){
	parse(genoMap, qualityMap);

	//recalibrate
	seqErrorModels.recalibrate(_bases);
	_sequenceAndQualitiesChanged = seqErrorModels.recalibrationChangesQualities();
};

void TAlignment::_parseBasesQualities(const TGenotypeMap & genoMap, const TQualityMap & qualityMap){
	//initialize
	_bases.resize(_cigar.lengthSequenced());
	_alignedPosition.resize(_cigar.lengthSequenced());
	int d = 0; //index regarding data structures and inside read
	int p = 0; //index regarding reference position (!= k for indels)

	//loop over cigar operations
	for(auto& cigarIter : _cigar){
		switch ( cigarIter.type ) {

			// for 'M', '=' or 'X': just copy
			case ('M') :
			case ('=') :
			case ('X') :
				//soft-clipped bases on 5' are before bamAlignment.Position
				for(unsigned int i=0; i<cigarIter.length; ++i, ++d, ++p){
					_bases[d].base = genoMap.toBase(_sequence[d]);
					_bases[d].originalQuality_phredInt = qualityMap.qualityToPhredInt(_qualities[d]);
					_bases[d].setAligned(true);
					_alignedPosition[d] = p;
				}
				break;

			//for 'S' - soft clip: ignore bases, but increase k
			case ('S') :
				//add bases to softclipped entries
				for(unsigned int i=0; i<cigarIter.length; ++i, ++d){
					//soft-clipped bases on 5' are before bamAlignment.Position
					//need to initialize quality for quality filter and bases for context
					_bases[d].base = genoMap.toBase(_sequence[d]);
					_bases[d].originalQuality_phredInt = qualityMap.qualityToPhredInt(_qualities[d]);
					_bases[d].setAligned(false);
					_alignedPosition[d] = -1;
				}
				break;

			//for 'I' - insertion: copy bases, but put aligned pos to -1
			case ('I')      :
				for(unsigned int i=0; i<cigarIter.length; ++i, ++d){
					_bases[d].base = genoMap.toBase(_sequence[d]);
					_bases[d].originalQuality_phredInt = qualityMap.qualityToPhredInt(_qualities[d]);
					_bases[d].setAligned(false);
					_alignedPosition[d] = -1;
				}
				break;


			// for 'D' - deletion: just add to position
			case ('D') :
				p += cigarIter.length;
				break;

			// for 'N' - skipped region in reference: only advance reference position
			case ('N') :
				p += cigarIter.length;
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
	_lastAlignedPositionWithRespectToRef = *this + (p - 1);
	_lastAlignedPos = p - 1; //why -1? -> same reason as above
};

void TAlignment::_setDistancesFromEnds(){
	//Set distances in ORIGINAL FRAGMENT (i.e. 5' end is where sequencing started, NOT how it aligns to reference)
	int length = _cigar.lengthSequenced();

	//is it paired-end?
	if(_flags.isProperPair()){
		if(_flags.isReverseStrand()){
			//reverse (can be either first or second mate, but it's the one that comes second in bam file)
			//and distance from 5' is given as f(end of fragment) = f(len - pos - 1)
			//hence distance from 3' is given by f(dist since beginning of fragment) = f(insert - len + pos)
			int k = abs(_fragmentLength) - (_cigar.lengthSequenced() - _cigar.lengthSoftClippedRight());
			int l = _cigar.lengthSequenced() - 1 - _cigar.lengthSoftClippedRight();
			for(int pos=0; pos<length; ++pos){
				_bases[pos].distFrom5Prime = l - pos; //dist from 5'
				_bases[pos].distFrom3Prime = k + pos; //dist from 3'
			}
		} else {
			//forward (can be either first or second mate, but it's the one that comes first in bam file)
			//Hence distance from 5' is given as a function of pos
			//And distance from 3' is given by (length of fragment) - pos -1
			//NOTE! we ignore indels when calculating distance from 5' since we can not know this info.
			//Luckily, this has only minimal effect since these distances are far from fragment ends
			for(int pos=0; pos<length; ++pos){
				_bases[pos].distFrom5Prime = pos - _cigar.lengthSoftClippedLeft(); //dist from 5'
				_bases[pos].distFrom3Prime = _fragmentLength - 1 - pos + _cigar.lengthSoftClippedLeft(); //dist from 3'
			}
		}
	} else {
		//treat as single end
		int l = length - 1;
		if(_flags.isReverseStrand()){
			//not in pair & reverse
			//Hence distance from 3' is just pos
			//And distance from 5' is just len - pos - 1
			for(int pos=0; pos<length; ++pos){
				_bases[pos].distFrom5Prime = l - pos - _cigar.lengthSoftClippedRight(); //dist from 5'
				_bases[pos].distFrom3Prime = pos - _cigar.lengthSoftClippedLeft(); //dist from 3'
			}
		} else {
			//not in pair & forward
			//Hence distance from 5' is just pos
			//And distance from 3' is given by len - pos - 1
			for(int pos=0; pos<length; ++pos){
				_bases[pos].distFrom5Prime = pos - _cigar.lengthSoftClippedLeft(); //dist from 5'
				_bases[pos].distFrom3Prime = l - pos - _cigar.lengthSoftClippedRight(); //dist from 3'
			}
		}
	}
};

void TAlignment::_fillContext(){
	if(_flags.isReverseStrand()){
		//reverse
		for(int d=0; d<(_cigar.lengthSequenced()-1); ++d){
			_bases[d].context = _bases[d+1].base;
		}
		_bases[_cigar.lengthSequenced()-1].context = N;
	} else {
		//forward
		_bases[0].context = N;
		for(int d=1; d<_cigar.lengthSequenced(); ++d)
			_bases[d].context = _bases[d-1].base;
	}
};

void TAlignment::addReference(TFastaBuffer & fasta){
	fasta.fill(*this, _lastAlignedPositionWithRespectToRef, _referenceSequence);
	_hasReference = true;
};

void TAlignment::setSequenceQualities(const TCigar Cigar, const std::string Sequence, const std::string Qualities){
	if(Cigar.lengthSequenced() != Sequence.length() || Cigar.lengthSequenced() != Qualities.length()){
		throw "Failed to set sequence and qualities of TAlignment: length of CIGAR, Sequences and Qualities does not match!";
	}
	_cigar = Cigar;
	_sequence = Sequence;
	_qualities = Qualities;
};

void TAlignment::setReadGroup(const uint16_t readGroupId){
	_readGroupID = readGroupId;
	_changed = true;
};

//--------------------------------------
//getters
//--------------------------------------
bool TAlignment::isAlignedAtInternalPos(const uint32_t internalPosition) const{
	return _alignedPosition[internalPosition] >= 0;
};

char TAlignment::referenceAtInternalPos(const uint32_t internalPosition) const{
	//Note: does not check if reference exists!
	return _referenceSequence[_alignedPosition[internalPosition]];
};

uint32_t TAlignment::positionInRef(const uint32_t internalPosition) const{
	//only makes sense if position is aligned!
	return position() + _alignedPosition[internalPosition];
};

uint16_t TAlignment::parsedLength() const{
	if(_parsed){
		return _cigar.lengthSequenced();
	} else {
		return 0;
	}
};

void TAlignment::_updateSequenceAndQualities(const TGenotypeMap & genoMap, const TQualityMap & qualMap) const{
	if(_sequenceAndQualitiesChanged){
		//update according to what is stored in bases
		_sequence.clear();
		_qualities.clear();

		for(auto& b : _bases){
			_sequence += genoMap.baseToChar[b.base];
			_qualities += (char) qualMap.phredIntToQuality(b.recalibratedQualityAsPhredInt);
		}

		_sequenceAndQualitiesChanged = false;
	}
};

std::string TAlignment::sequence(const TGenotypeMap & genoMap, const TQualityMap & qualMap) const{
	_updateSequenceAndQualities(genoMap, qualMap);
	return _sequence;
};

std::string TAlignment::qualities(const TGenotypeMap & genoMap, const TQualityMap & qualMap) const{
	_updateSequenceAndQualities(genoMap, qualMap);
	return _qualities;
};

//--------------------------------------------
//filters and other functions to modify data
//--------------------------------------------
void TAlignment::filterForBaseQuality(TQualityFilter & qualFilter){
	//set base to N if outside quality filter
	for(auto& b : _bases){
		if(!qualFilter.pass(b.recalibratedQualityAsPhredInt)){
			b.base = N;
			b.recalibratedQualityAsPhredInt = 0;
		}
	}

	_sequenceAndQualitiesChanged = true;
	_changed = true;
};

void TAlignment::filterForContext(const std::map<BaseContext,int> & ignoreTheseContexts, const TGenotypeMap & genoMap){
	//set base to N if outside quality filter
	for(auto& b : _bases){
		BaseContext c = genoMap.toContext(b.base, b.context);
		if(ignoreTheseContexts.find(c) != ignoreTheseContexts.end()){
			b.base = N;
			b.recalibratedQualityAsPhredInt = 0;
		}
	}

	_sequenceAndQualitiesChanged = true;
	_changed = true;
};

void TAlignment::trimRead(const int & trimmingLength3Prime, const int & trimmingLength5Prime){
	for(auto& b : _bases){
		if(b.distFrom3Prime < trimmingLength3Prime || b.distFrom5Prime < trimmingLength5Prime){
			b.base = N;
			b.recalibratedQualityAsPhredInt = 0;
		}
	}

	_sequenceAndQualitiesChanged = true;
	_changed = true;
};

void TAlignment::removeSoftClippedBases(){
	//check if there is softclipping
	if(_cigar.lengthSoftClipped() > 0){
		auto bIter = _bases.begin();
		for(auto& cigarIter : _cigar){
			if(cigarIter.type == 'S'){
				//remove bases
				bIter = _bases.erase(bIter, bIter + cigarIter.length);
			} else {
				//just advance position
				bIter += cigarIter.length;
			}
		}

		//update cigar and length
		_cigar.removeSoftClips();

		//set has changed
		_sequenceAndQualitiesChanged = true;
		_changed = true;
	}
};

void TAlignment::binQualityScores(TQualityMap & qualityMap){
	//make sure read is parsed
	if(!_parsed) throw "Read was not parsed!";

	//bin quality scores as done by Illumina
	for(auto& b : _bases){
		b.recalibratedQualityAsPhredInt = qualityMap.phredIntToIlluminaPhredInt(b.recalibratedQualityAsPhredInt);
	}

	_sequenceAndQualitiesChanged = true;
	_changed = true;
};

void TAlignment::recalibrateWithPMD(const GenotypeLikelihoods::TGenotypeLikelihoodCalculator & GLCalculator){
	GLCalculator.recalibrateWithPMD(_bases);
	_sequenceAndQualitiesChanged = true;
	_changed = true;
};

void TAlignment::setIsProperPair(const bool & ok){
	_flags.setIsProperPair(ok);
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

void TAlignment::downsampleAlignment(const double fractionToKeep, TRandomGenerator& randomGenerator){
	for(auto& b : _bases){
		double r = randomGenerator.getRand();
		if(r > fractionToKeep){
			b.base = N;
			b.recalibratedQualityAsPhredInt = 0;
		}
	}
	_sequenceAndQualitiesChanged = true;
	_changed = true;
};

//--------------------------------------------
//functions to write / print alignment
//--------------------------------------------
void TAlignment::print(TGenotypeMap & genoMap, TQualityMap & qualMap){
	std::cout << "NAME:\t" << _name << std::endl;
	std::cout << "LEN:\t" << _bases.size() << std::endl;

	//print bases
	std::cout << "SEQ:\t";
	for(auto& b : _bases){
		std::cout << genoMap.getBaseAsChar(b.base);
	}
	std::cout << std::endl;

	//print qualities
	std::cout << "QUAL:\t";
	for(auto& b : _bases){
		std::cout << (char) qualMap.phredIntToQuality(b.recalibratedQualityAsPhredInt);
	std::cout << std::endl;

	//print aligned pos
	std::cout << "POS:\t";
	for(size_t d=0; d<_bases.size(); ++d){
		if(d>0) std::cout << ",";
		if(_bases[d].isAligned())
			std::cout << _alignedPosition[d];
		else
			std::cout << "-";
	}
	std::cout << std::endl;

	//print dist from 3'
	std::cout << "dist 3':\t";
	bool first = true;
	for(auto& b : _bases){
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
	for(auto& b : _bases){
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
