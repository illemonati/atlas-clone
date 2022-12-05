/*
 * TAlignment.cpp
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */

#include "TAlignment.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "TBamFilter.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/strongTypes.h"
#include <algorithm>
#include <iterator>
#include <iostream>
#include <math.h>
#include <memory>
#include <stdexcept>
#include <stdlib.h>

namespace BAM {

void TAlignment::clear() {
	TGenomePosition::clear();
	_name.clear();
	_flags.reset();
	_mappingQuality = 0;
	_cigar.clear();
	_mateGenomicPosition.clear();
	_insertSize_TLEN = 0;
	_readGroupID    = 0;
	_fragmentLength = 0;

	_lastAlignedPos = 0;

	// booleans
	_empty   = true;
	_parsed  = false;

	// sequence and qualities. Mutable so that they can be recreated from bases even for const TAlignment
	_sequence.clear();
	_qualities.clear();
	_sequenceAndQualitiesChanged = false;

	// per base data
	_bases.clear();
	_alignedPosition.clear();

	// reference
	_referenceSequence.clear();
}

//--------------------------------------
// functions to fill alignment
//--------------------------------------
// function used by TBamFile to fill alignment
void TAlignment::fill(const std::string &Name, const TSamFlags &Flags, uint32_t RefID, uint32_t Position,
					  uint16_t MappingQuality, const TCigar &Cigar, uint32_t MateRefID, uint32_t MatePosition,
					  const int32_t &InsertSize_TLEN, const std::string &Sequence, const std::string &Qualities,
					  uint16_t ReadGroupId) {

	// empty alignment
	_lastAlignedPos = 0;
	_parsed  = false;
	_sequenceAndQualitiesChanged = false;
	_bases.clear();
	_alignedPosition.clear();
	_referenceSequence.clear();

	// copy data
	_name  = Name;
	_flags = Flags;
	move(RefID, Position);
	_mappingQuality  = MappingQuality;
	_insertSize_TLEN = InsertSize_TLEN;
	_sequence        = Sequence;
	_qualities       = Qualities;
	_readGroupID     = ReadGroupId;
	_empty           = false;

	if (_flags.isPaired()) {
		_mateGenomicPosition.move(MateRefID, MatePosition);
	} else {
		_mateGenomicPosition.move(0, 0); // 0 is not paired
	}

	_setCigar(Cigar);
};

void TAlignment::_setCigar(const TCigar &Cigar){
	_cigar = Cigar;
	_parsed = false;
	// set fragment length
	if (_flags.isProperPair()) {
		_fragmentLength = abs(_insertSize_TLEN) + _cigar.lengthInserted() - _cigar.lengthDeleted();
	} else {
		_fragmentLength = _cigar.lengthSequenced();
	}
};

void TAlignment::_parseBasesQualities() {
	using genometools::char2base;
	if (_sequence.size() != _qualities.size()) {
		throw std::runtime_error(
			"void TAlignment::_parseBasesQualities(): Sequence and Qualities are of different legth!");
	}
	// initialize
	_bases.resize(_cigar.lengthRead());
	//_alignedPosition.resize(_cigar.lengthRead());
	_alignedPosition.clear();
	_alignedPosition.reserve(_cigar.lengthRead());
	int d = 0; // index regarding data structures and inside read
	int p = 0; // index regarding reference position (!= d for soft clipping & indels)

	// loop over cigar operations
	for (auto &cigarIter : _cigar) {
		switch (cigarIter.type) {

		// for 'M', '=' or 'X': just copy
		case ('M'):
		case ('='):
		case ('X'):
			// soft-clipped bases on 5' are before bamAlignment.Position
			for (unsigned int i = 0; i < cigarIter.length; ++i, ++d, ++p) {
				_bases[d].base                     = char2base(_sequence[d]);
				_bases[d].originalQuality_phredInt = genometools::BaseQuality(_qualities[d]);
				_bases[d].setAligned(true);
				_alignedPosition.push_back(p);
			}
			_lastAlignedPos = d - 1; // Note: for loop ends with d one too large
			break;

		// for 'S' - soft clip: copy by set aligned = false
		case ('S'):
			// add bases to softclipped entries
			for (unsigned int i = 0; i < cigarIter.length; ++i, ++d) {
				// soft-clipped bases on 5' are before bamAlignment.Position
				// need to initialize quality for quality filter and bases for context
				_bases[d].base                     = char2base(_sequence[d]);
				_bases[d].originalQuality_phredInt = genometools::BaseQuality(_qualities[d]);
				_bases[d].setAligned(false);
				_alignedPosition.push_back(-1);
			}
			break;

		// for 'I' - insertion: copy bases, but put aligned  = false
		case ('I'):
			for (unsigned int i = 0; i < cigarIter.length; ++i, ++d) {
				_bases[d].base                     = char2base(_sequence[d]);
				_bases[d].originalQuality_phredInt = genometools::BaseQuality(_qualities[d]);
				_bases[d].setAligned(false);
				_alignedPosition.push_back(-1);
			}
			break;

		// for 'D' - deletion: just add to position
		case ('D'): p += cigarIter.length; break;

		// for 'N' - skipped region in reference: only advance reference position
		case ('N'): p += cigarIter.length; break;

		// for 'H' or 'P' - hard clip: do nothing as these bases are not present in SEQ
		case ('H'):
		case ('P'): break;

		// invalid CIGAR op-code
		default: throw(std::string) "CIGAR operation '" + cigarIter.type + "' not supported!";
		}
	}

	// update length and last aligned position
	_refSize = p;

	// then update distances from ends
	_setDistancesFromEnds();

	// fill context for each base
	_fillContext();

	// set mapping quality and whether read is first or second
	for (auto &b : _bases) {
		b.readGroupID    = _readGroupID;
		b.mappingQuality = _mappingQuality;
		b.fragmentLength = _fragmentLength;
		b.setSecondMate(_flags.isSecondMate());
		b.setReverseStrand(_flags.isReverseStrand());
	}

	_parsed                      = true;
	_sequenceAndQualitiesChanged = false;
};

void TAlignment::_setQualitiesNoRecal() {
	// No recal: set recalibrated quality = original quality
	for (auto &b : _bases) { b.recalibratedQualityAsPhredInt = b.originalQuality_phredInt; }
};

void TAlignment::_setDistancesFromEnds() {
	// Set distances in ORIGINAL FRAGMENT (i.e. 5' end is where sequencing started, NOT how it aligns to reference)
	const int length = _cigar.lengthSequenced();

	// is it paired-end?
	if (_flags.isProperPair()) {
		if (_flags.isReverseStrand()) {
			// reverse (can be either first or second mate, but it's the one that comes second in bam file)
			// and distance from 5' is given as f(end of fragment) = f(len - pos - 1)
			// hence distance from 3' is given by f(dist since beginning of fragment) = f(insert - len + pos)
			const int k = _fragmentLength - (_cigar.lengthSequenced() - _cigar.lengthSoftClippedRight());
			const int l = _cigar.lengthSequenced() - 1 - _cigar.lengthSoftClippedRight();
			for (int pos = 0; pos < length; ++pos) {
				_bases[pos].distFrom5Prime = l - pos; // dist from 5'
				_bases[pos].distFrom3Prime = k + pos; // dist from 3'
			}
		} else {
			// forward (can be either first or second mate, but it's the one that comes first in bam file)
			// Hence distance from 5' is given as a function of pos
			// And distance from 3' is given by (length of fragment) - pos -1
			// NOTE! we ignore indels when calculating distance from 5' since we can not know this info.
			// Luckily, this has only minimal effect since these distances are far from fragment ends
			for (int pos = 0; pos < length; ++pos) {
				_bases[pos].distFrom5Prime = pos - _cigar.lengthSoftClippedLeft();                       // dist from 5'
				_bases[pos].distFrom3Prime = _fragmentLength - 1 - pos + _cigar.lengthSoftClippedLeft(); // dist from 3'
			}
		}
	} else {
		// treat as single end
		const int l = length - 1;
		if (_flags.isReverseStrand()) {
			// not in pair & reverse
			// Hence distance from 3' is just pos
			// And distance from 5' is just len - pos - 1
			for (int pos = 0; pos < length; ++pos) {
				_bases[pos].distFrom5Prime = l - pos - _cigar.lengthSoftClippedRight(); // dist from 5'
				_bases[pos].distFrom3Prime = pos - _cigar.lengthSoftClippedLeft();      // dist from 3'
			}
		} else {
			// not in pair & forward
			// Hence distance from 5' is just pos
			// And distance from 3' is given by len - pos - 1
			for (int pos = 0; pos < length; ++pos) {
				_bases[pos].distFrom5Prime = pos - _cigar.lengthSoftClippedLeft();      // dist from 5'
				_bases[pos].distFrom3Prime = l - pos - _cigar.lengthSoftClippedRight(); // dist from 3'
			}
		}
	}
};

void TAlignment::_fillContext() {
	using namespace genometools;
	if (_flags.isReverseStrand()) {
		// reverse
		for (size_t d = 0; d < (_cigar.lengthSequenced() - 1); ++d) {
			_bases[d].previousBase = _bases[d + 1].base;
		}
		_bases[_cigar.lengthSequenced() - 1].previousBase = Base::N;
	} else {
		// forward
		_bases[0].previousBase = Base::N;
		for (size_t d = 1; d < _cigar.lengthSequenced(); ++d)
			_bases[d].previousBase = _bases[d - 1].base;
	}
};

void TAlignment::parse() {
	_parseBasesQualities();
	_setQualitiesNoRecal();
};

void TAlignment::parse(const GenotypeLikelihoods::SequencingError::TModels &seqErrorModels) {
	_parseBasesQualities();

	// recalibrate
	seqErrorModels.recalibrate(_bases);
	_sequenceAndQualitiesChanged = seqErrorModels.recalibrationChangesQualities();
};

void TAlignment::addReference(const genometools::TFastaReader &fasta) {
	const auto window = genometools::TGenomeWindow(*this, _lastAlignedPositionWithRespectToRef);
	const auto view = fasta.view(window);
	_referenceSequence.clear();
	std::copy(view.begin(), view.end(), std::back_inserter(_referenceSequence));
};

void TAlignment::setSequenceQualities(const TCigar &Cigar, const std::vector<genometools::Base> &Sequence,
									  const std::vector<genometools::PhredIntProbability> &Qualities) {
	if (Cigar.lengthRead() != Sequence.size() || Cigar.lengthRead() != Qualities.size()) {
		throw std::runtime_error(
			"void TAlignment::setSequenceQualities(const TCigar & Cigar, const std::vector<Base> & Sequence, const "
			"std::vector<PhredIntProbability> & Qualities): length of CIGAR, Sequences and Qualities do not match!");
	}
	_setCigar(Cigar);

	// parse bases and qualities
	_sequence.clear();
	_sequence.reserve(Sequence.size());
	_qualities.clear();
	_qualities.reserve(Sequence.size());
	for (size_t i = 0; i < Sequence.size(); ++i) {
		_sequence.push_back(genometools::base2char(Sequence[i]));
		_qualities.push_back(static_cast<char>(genometools::BaseQuality(Qualities[i])));
	}
	_parseBasesQualities();	
	_setQualitiesNoRecal();
	_sequenceAndQualitiesChanged = true; // will trigger that the strings are read form the bases
};

void TAlignment::setReadGroup(const uint16_t readGroupId) {
	_readGroupID = readGroupId;
};
//get overlap length from previous function and pass it on to cigar constructor
void TAlignment::merge(uint16_t overlapLength, size_t &mappedBasesClipped) {
	_setCigar(TCigar(cigar(), overlapLength, !isReverseStrand(), mappedBasesClipped));
}

//--------------------------------------
// getters
//--------------------------------------
bool TAlignment::isAlignedAtInternalPos(const uint32_t internalPosition) const {
	return _alignedPosition[internalPosition] >= 0;
};

uint32_t TAlignment::getLastInternalPos() const {
	return (_alignedPosition.size()-1);
}

genometools::Base TAlignment::referenceAtInternalPos(uint32_t internalPosition) const {
 	assert(isAlignedAtInternalPos(internalPosition));
	assert((size_t)_alignedPosition[internalPosition] < _referenceSequence.size());
 	return _referenceSequence[_alignedPosition[internalPosition]];
};

genometools::TGenomePosition TAlignment::positionInRef(uint32_t internalPosition) const {
	// only makes sense if position is aligned!
	return *this + _alignedPosition[internalPosition];
};

void TAlignment::_updateSequenceAndQualities() const {
	if (_sequenceAndQualitiesChanged) {
		// update according to what is stored in bases
		_sequence.resize(_bases.size());
		_qualities.resize(_bases.size());

		for (size_t b = 0; b < _bases.size(); ++b) {
			_sequence[b]  = genometools::base2char(_bases[b].base);
			_qualities[b] = static_cast<char>(genometools::BaseQuality(_bases[b].recalibratedQualityAsPhredInt));
		}

		_sequenceAndQualitiesChanged = false;
	}
};

std::string TAlignment::sequence() const {
	_updateSequenceAndQualities();
	return _sequence;
};

std::string TAlignment::qualities() const {
	_updateSequenceAndQualities();
	return _qualities;
};

//--------------------------------------------
// filters and other functions to modify data
//--------------------------------------------
void TAlignment::filter(const TBaseFilter &Filter) {
	if (Filter) {
		// set quality = 0 and base = N if outside quality filter
		for (auto &b : _bases) {
			if (!Filter.pass(b)) {
				b.base                          = genometools::Base::N;
				b.recalibratedQualityAsPhredInt = 0;
			}
		}
	}
};

void TAlignment::trimRead(int trimmingLength3Prime, int trimmingLength5Prime) {
	for (auto &b : _bases) {
		if (b.distFrom3Prime < trimmingLength3Prime || b.distFrom5Prime < trimmingLength5Prime) {
			b.base                          = genometools::Base::N;
			b.recalibratedQualityAsPhredInt = 0;
		}
	}

	_sequenceAndQualitiesChanged = true;
};

void TAlignment::removeSoftClippedBases() {
	// make sure read is parsed
	if (!_parsed) throw std::runtime_error("void TAlignment::removeSoftClippedBases(): Read was not parsed!");

	// check if there is softclipping
	if (_cigar.lengthSoftClipped() > 0) {
		auto bIter = _bases.begin();
		for (auto &cigarIter : _cigar) {
			if (cigarIter.type == 'S') {
				// remove bases
				bIter = _bases.erase(bIter, bIter + cigarIter.length);
			} else if (cigarIter.type == 'M' || cigarIter.type == '=' || cigarIter.type == 'X' ||
					   cigarIter.type == 'I') {
				// just advance position
				bIter += cigarIter.length;
			}
		}

		// update cigar and length
		_cigar.removeSoftClips();

		// set has changed
		_sequenceAndQualitiesChanged = true;
	}
};

void TAlignment::binQualityScoresIllumina() {
	// make sure read is parsed
	if (!_parsed)
		throw std::runtime_error("void TAlignment::binQualityScores(TQualityMap & qualityMap): Read was not parsed!");

	// bin quality scores as done by Illumina
	for (auto &b : _bases) { b.recalibratedQualityAsPhredInt.makeIllumina(); }

	_sequenceAndQualitiesChanged = true;
};

void TAlignment::recalibrateWithPMD(const GenotypeLikelihoods::TGenotypeLikelihoodCalculator &GLCalculator) {
	GLCalculator.recalibrateWithPMD(_bases);
	_sequenceAndQualitiesChanged = true;
};

void TAlignment::setIsProperPair(const bool &ok) { _flags.setIsProperPair(ok); };

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

void TAlignment::downsampleAlignment(const coretools::Probability &fractionToKeep,
									 coretools::TRandomGenerator &randomGenerator) {
	for (auto &b : _bases) {
		double r = randomGenerator.getRand();
		if (r > fractionToKeep) {
			b.base                          = genometools::Base::N;
			b.recalibratedQualityAsPhredInt = 0;
		}
	}
	_sequenceAndQualitiesChanged = true;
};

//--------------------------------------------
// functions to write / print alignment
//--------------------------------------------
void TAlignment::print() const {
	std::cout << std::endl << "NAME:\t" << _name << std::endl;
	std::cout << "LEN:\t" << _bases.size() << std::endl;

	// print bases
	std::cout << "SEQ:\t";
	for (auto &b : _bases) { std::cout << b.base; }
	std::cout << std::endl;

	// print qualities
	std::cout << "QUAL:\t";
	for (auto &b : _bases) { std::cout << genometools::BaseQuality(b.recalibratedQualityAsPhredInt); }
	std::cout << std::endl;

	// print aligned pos
	std::cout << "POS:\t";
	for (size_t d = 0; d < _bases.size(); ++d) {
		if (d > 0) std::cout << ",";
		if (_bases[d].isAligned())
			std::cout << _alignedPosition[d];
		else
			std::cout << "-";
	}
	std::cout << std::endl;

	// print dist from 3'
	std::cout << "dist 3':\t";
	bool first = true;
	for (auto &b : _bases) {
		if (first) {
			first = false;
		} else {
			std::cout << ",";
		}
		std::cout << b.distFrom3Prime;
	}
	std::cout << std::endl;

	// print dist from 5'
	std::cout << "dist 5':\t";
	first = true;
	for (auto &b : _bases) {
		if (first) {
			first = false;
		} else {
			std::cout << ",";
		}
		std::cout << b.distFrom5Prime;
	}
	std::cout << std::endl;
};

}; // namespace BAM
