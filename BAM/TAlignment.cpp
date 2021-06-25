/*
 * TAlignment.cpp
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */

#include "TAlignment.h"

namespace BAM{

TAlignment::TAlignment(const uint32_t& RefID, const uint32_t& Position):TGenomePosition(RefID, Position){
	_initialize();
};

TAlignment::TAlignment(const TGenomePosition & other):TGenomePosition(other){
	_initialize();
};

TAlignment::TAlignment(){
	_initialize();
};

void TAlignment::_initialize(){
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
void TAlignment::fill(const	std::string & Name,
		  const TSamFlags & Flags,
		  const uint32_t & RefID,
		  const uint32_t & Position,
		  const uint16_t & MappingQuality,
		  const TCigar & Cigar,
		  const uint32_t & MateRefID,
		  const uint32_t & MatePosition,
		  const int32_t & InsertSize_TLEN,
		  const std::string & Sequence,
		  const std::string & Qualities,
		  const uint16_t & ReadGroupId){

	//empty alignment
	clear();

	//copy data
	_name = Name;
	_flags = Flags;
	move(RefID, Position);
	_mappingQuality = MappingQuality;
	_cigar = Cigar;
	_insertSize_TLEN = InsertSize_TLEN;
	_sequence = Sequence;
	_qualities = Qualities;
	_readGroupID = ReadGroupId;
	_empty = false;

	if (_flags.isPaired()){
        _mateGenomicPosition.move(MateRefID, MatePosition);
    } else {
        _mateGenomicPosition.move(0, 0); // 0 is not paired
    }
	//set fragment length
	if(_flags.isProperPair()){
		_fragmentLength = abs(_insertSize_TLEN) + _cigar.lengthInserted() - _cigar.lengthDeleted();
	} else {
		_fragmentLength = _cigar.lengthSequenced();
	}
};

void TAlignment::parse(){
	//first parse bases and qualities
	_parseBasesQualities(_sequence, _qualities);
};

void TAlignment::parse(const GenotypeLikelihoods::TSequencingErrorModels & seqErrorModels){
	parse();

	//recalibrate
	seqErrorModels.recalibrate(_bases);
	_sequenceAndQualitiesChanged = seqErrorModels.recalibrationChangesQualities();
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
			_bases[d].context.set(_bases[d+1].base, _bases[d].base);
		}
		_bases[_cigar.lengthSequenced()-1].context.set(genometools::N, _bases[_cigar.lengthSequenced()-1].base);
	} else {
		//forward
		_bases[0].context.set(genometools::N, _bases[0].base);
		for(int d=1; d<_cigar.lengthSequenced(); ++d)
			_bases[d].context.set(_bases[d-1].base, _bases[d].base);
	}
};

void TAlignment::addReference(TFastaBuffer & fasta){
	fasta.fill(*this, _lastAlignedPositionWithRespectToRef, _referenceSequence);
	_hasReference = true;
};

void TAlignment::setSequenceQualities(const TCigar & Cigar, const std::vector<genometools::Base> & Sequence, const std::vector<genometools::PhredIntProbability> & Qualities){
	if(Cigar.lengthRead() != Sequence.size() || Cigar.lengthRead() != Qualities.size()){
		throw std::runtime_error("void TAlignment::setSequenceQualities(const TCigar & Cigar, const std::vector<Base> & Sequence, const std::vector<PhredIntProbability> & Qualities): length of CIGAR, Sequences and Qualities do not match!");
	}
	_cigar = Cigar;

	//parse bases and qualities
	_parseBasesQualities(_sequence, _qualities);

	_changed = true;
	_sequenceAndQualitiesChanged = true; //will trigger that the strings are read form the bases
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

BAM::Base TAlignment::referenceAtInternalPos(const uint32_t & internalPosition) const{
	if(!_hasReference){
		throw std::runtime_error("BAM::Base TAlignment::referenceAtInternalPos(const uint32_t internalPosition) const: alignment has no reference!");
	}
	return _referenceSequence[_alignedPosition[internalPosition]];
};

TGenomePosition TAlignment::positionInRef(const uint32_t & internalPosition) const{
	//only makes sense if position is aligned!
	return *this + _alignedPosition[internalPosition];
};

uint16_t TAlignment::parsedLength() const{
	if(_parsed){
		return _cigar.lengthRead();
	} else {
		return 0;
	}
};

void TAlignment::_updateSequenceAndQualities() const{
	if(_sequenceAndQualitiesChanged){
		//update according to what is stored in bases
		_sequence.resize(_bases.size());
		_qualities.resize(_bases.size());

		for(auto b=0; b < _bases.size(); ++b){
			_sequence[b] = (char) _bases[b].base;
			_qualities[b] = (char) genometools::BaseQuality(_bases[b].recalibratedQualityAsPhredInt);
		}

		_sequenceAndQualitiesChanged = false;
	}
};

std::string TAlignment::sequence() const{
	_updateSequenceAndQualities();
	return _sequence;
};

std::string TAlignment::qualities() const{
	_updateSequenceAndQualities();
	return _qualities;
};

//--------------------------------------------
//filters and other functions to modify data
//--------------------------------------------
void TAlignment::filter(const TBaseFilter & Filter){
	if(Filter){
		//set quality = 0 and base = N if outside quality filter
		for(auto& b : _bases){
			if(!Filter.pass(b)){
				b.base = genometools::N;
				b.recalibratedQualityAsPhredInt = 0;
			}
		}
	}
};

void TAlignment::trimRead(const int & trimmingLength3Prime, const int & trimmingLength5Prime){
	for(auto& b : _bases){
		if(b.distFrom3Prime < trimmingLength3Prime || b.distFrom5Prime < trimmingLength5Prime){
			b.base = genometools::N;
			b.recalibratedQualityAsPhredInt = 0;
		}
	}

	_sequenceAndQualitiesChanged = true;
	_changed = true;
};

void TAlignment::removeSoftClippedBases(){
	//make sure read is parsed
	if(!_parsed) throw std::runtime_error("void TAlignment::removeSoftClippedBases(): Read was not parsed!");

	//check if there is softclipping
	if(_cigar.lengthSoftClipped() > 0){
		auto bIter = _bases.begin();
		for(auto& cigarIter : _cigar){
			if(cigarIter.type == 'S'){
				//remove bases
				bIter = _bases.erase(bIter, bIter + cigarIter.length);
			} else if(cigarIter.type == 'M' || cigarIter.type == '=' || cigarIter.type == 'X' || cigarIter.type == 'I'){
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

void TAlignment::binQualityScoresIllumina(){
	//make sure read is parsed
	if(!_parsed) throw std::runtime_error("void TAlignment::binQualityScores(TQualityMap & qualityMap): Read was not parsed!");

	//bin quality scores as done by Illumina
	for(auto& b : _bases){
		b.recalibratedQualityAsPhredInt.makeIllumina();
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

void TAlignment::downsampleAlignment(const coretools::Probability & fractionToKeep, coretools::TRandomGenerator& randomGenerator){
	for(auto& b : _bases){
		double r = randomGenerator.getRand();
		if(r > fractionToKeep){
			b.base = genometools::N;
			b.recalibratedQualityAsPhredInt = 0;
		}
	}
	_sequenceAndQualitiesChanged = true;
	_changed = true;
};

//--------------------------------------------
//functions to write / print alignment
//--------------------------------------------
void TAlignment::print(){
	std::cout << std::endl << "NAME:\t" << _name << std::endl;
	std::cout << "LEN:\t" << _bases.size() << std::endl;

	//print bases
	std::cout << "SEQ:\t";
	for(auto& b : _bases){
		std::cout << b.base;
	}
	std::cout << std::endl;

	//print qualities
	std::cout << "QUAL:\t";
	for(auto& b : _bases){
		std::cout << genometools::BaseQuality(b.recalibratedQualityAsPhredInt);
	}
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
	first = true;
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
