/*
 * TAlignment.cpp
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */
#include "TAlignment.h"

#include "SequencingError/TModels.h"
#include "TErrorModels.h"
#include "TSequencedBase.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/TLogInt.h"
#include "coretools/Types/TPseudoInt.h"
#include "coretools/Types/probability.h"
#include "genometools/TFastaReader.h"
#include <cstdint>

namespace BAM {

using coretools::PhredInt;
using coretools::fromChar;

PhredInt makeIllumina(PhredInt value) noexcept {
	if (value < 35)
		return PhredInt(33);
	if (value < 43)
		return PhredInt(39);
	if (value < 53)
		return PhredInt(48);
	if (value < 58)
		return PhredInt(55);
	if (value < 63)
		return PhredInt(60);
	if (value < 68)
		return PhredInt(66);
	if (value < 72)
		return PhredInt(70);

	return PhredInt(73);
}

char makeIllumina(char Quality) noexcept {
	const auto phr    = fromChar(Quality);
	const auto phrIll = makeIllumina(phr);
	return coretools::toChar(phrIll);
}

void TAlignment::clear() {
	TGenomePosition::clear();
	_name.clear();
	_flags.reset();
	_mappingQuality = PhredInt::highest();
	_cigar.clear();
	_mateGenomicPosition.clear();
	_insertSize_TLEN = 0;
	_readGroupID     = 0;
	_bamID           = 0;
	_fragmentLength  = 0;

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
					  size_t BamID, size_t ReadGroupId) {

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
	MappingQuality   = std::min<uint16_t>(MappingQuality, PhredInt::max());
	_mappingQuality  = PhredInt(MappingQuality);
	_insertSize_TLEN = InsertSize_TLEN;
	_sequence        = Sequence;
	_qualities       = Qualities;
	_bamID           = BamID;
	_readGroupID     = ReadGroupId;
	_empty           = false;

	if (_flags.isPaired()) {
		_mateGenomicPosition.move(MateRefID, MatePosition);
	} else {
		_mateGenomicPosition.move(0, 0); // 0 is not paired
	}

	_setCigar(Cigar);
}

void TAlignment::_setCigar(const TCigar &Cigar){
	_cigar = Cigar;
	_parsed = false;
	// set fragment length
	if (_flags.isProperPair()) {
		_fragmentLength = abs(_insertSize_TLEN) + _cigar.lengthInserted() - _cigar.lengthDeleted();
	} else {
		_fragmentLength = _cigar.lengthSequenced();
	}
}

void TAlignment::_parseBasesQualities() {
	using genometools::char2base;
	using coretools::TLogInt;
	if (_sequence.size() != _qualities.size()) {
		DEVERROR("Sequence and Qualities are of different legth!");
	}
	// initialize
	const auto common = [&](){
		// set mapping quality and whether read is first or second
		TSequencedBase b{};
		b.readGroupID    = _readGroupID;
		b.bamID          = _bamID;
		b.mappingQuality = _mappingQuality;
		b.fragmentLength = TLogInt::fromLinear(_fragmentLength);
		b.set<Flags::SecondMate>(_flags.isSecondMate());
		b.set<Flags::ReversedStrand>(_flags.isReverseStrand());
		return b;
	}();
	_bases.assign(_cigar.lengthRead(), common);
	_alignedPosition.resize(_cigar.lengthRead());
	int d = 0; // index regarding data structures and inside read
	int p = 0; // index regarding reference position (!= d for soft clipping & indels)

	// loop over cigar operations
	for (auto &cigarIter : _cigar) {
		switch (cigarIter.type) {

		// for 'M', '=' or 'X': just copy
		case ('M'):
		case ('='):
		case ('X'):
			// soft-clipped bases on left are before bamAlignment.Position
			for (unsigned int i = 0; i < cigarIter.length; ++i, ++d, ++p) {
				_bases[d].base            = char2base(_sequence[d]);
				_bases[d].originalQuality = fromChar(_qualities[d]);
				_bases[d].set<Flags::Aligned>(true);
				_bases[d].set<Flags::SoftClipped>(false);
				_alignedPosition[d] = p;
			}
			_lastAlignedPos = d - 1; // Note: for loop ends with d one too large
			break;

		// for 'S' - soft clip: copy by set aligned = false
		case ('S'):
			// add bases to softclipped entries
			for (unsigned int i = 0; i < cigarIter.length; ++i, ++d) {
				// soft-clipped bases on 5' are before bamAlignment.Position
				// need to initialize quality for quality filter and bases for context
				_bases[d].base            = char2base(_sequence[d]);
				_bases[d].originalQuality = fromChar(_qualities[d]);
				_bases[d].set<Flags::Aligned>(false);
				_bases[d].set<Flags::SoftClipped>(true);
				_alignedPosition[d] = -1;
			}
			break;

		// for 'I' - insertion: copy bases, but put aligned  = false
		case ('I'):
			for (unsigned int i = 0; i < cigarIter.length; ++i, ++d) {
				_bases[d].base            = char2base(_sequence[d]);
				_bases[d].originalQuality = fromChar(_qualities[d]);
				_bases[d].set<Flags::Aligned>(false);
				_bases[d].set<Flags::SoftClipped>(false);
				_alignedPosition[d] = -1;
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
		default: UERROR("CIGAR operation '", cigarIter.type, "' not supported!");
		}
	}

	// set mapping quality and whether read is first or second

	// update length and last aligned position
	_refSize = p;

	// then update distances from ends
	_setDistancesFromEnds();

	// fill context for each base
	_fillContext();


	_parsed                      = true;
	_sequenceAndQualitiesChanged = false;
}

void TAlignment::_setQualitiesNoRecal() {
	// No recal: set recalibrated quality = original quality
	for (auto &b : _bases) { b.recalQuality = b.originalQuality; }
}

void TAlignment::_setDistancesFromEnds() {
	using coretools::TPseudoInt;
	// Set distances in ORIGINAL FRAGMENT (i.e. 5' end is where sequencing started, NOT how it aligns to reference)
	const int length  = _cigar.lengthSequenced();
	const int l_m1_ps = length - 1 + _cigar.lengthSoftClippedLeft();

	// is it paired-end?
	if (_flags.isProperPair()) {
		if (_flags.isReverseStrand()) {
			// Paired-end reverse
			//    FFFFFFFFFF       : fragmentLenght   = 10
			//  3'  LL---x--RR  5' : readLength       = softclippedLeft + lengthSequenced + softClippedRight
			//  3'  ..---x--..  5' : lengthSequenced  = 6
			//  3'  LL...x....  5' : softClippledLeft = 2
			//  3'  .....x..RR  5' : softClippledRight= 2
			// pos: 0123456789
			// d5':   543210       : d5 = lengthSequenced - 1 - pos + softClippedLeft               = 6 - 1 - 5 + 2   = 2
			// d3':   456789       : d3 = pos + fragmentLenght - lengthSequenced  - softClippedLeft = 5 + 10 - 6 - 2  = 7
			const int f_ml_ms = _fragmentLength - _cigar.lengthSequenced() - _cigar.lengthSoftClippedLeft();
			for (size_t pos = _cigar.lengthSoftClippedLeft(); pos < length + _cigar.lengthSoftClippedLeft(); ++pos) {
				_bases[pos].distFrom5 = TPseudoInt::fromLinear(l_m1_ps - pos);
				_bases[pos].distFrom3 = TPseudoInt::fromLinear(pos + f_ml_ms);
			}
		} else {
			// Paired-end forward
			//        FFFFFFFFFF   : fragmentLenght   = 10
			//  5'  LL---x--RR  3' : readLength       = softclippedLeft + lengthSequenced + softClippedRight
			//  5'  ..---x--..  3' : lengthSequenced  = 6
			//  5'  LL...x....  3' : softClippledLeft = 2
			//  5'  .....x..RR  3' : softClippledRight= 2
			// pos: 0123456789
			// d5':   012345       : d5 = pos - softClippedLeft                       = 5 - 2          = 3
			// d3':   987654       : d3 = fragmentLength - 1 - pos + softClippedLeft  = 10 - 1 - 5 + 2 = 6
			const int f_m1_ps = _fragmentLength - 1 + _cigar.lengthSoftClippedLeft();
			for (size_t pos = _cigar.lengthSoftClippedLeft(); pos < length + _cigar.lengthSoftClippedLeft(); ++pos) {
				_bases[pos].distFrom5 = TPseudoInt::fromLinear(pos - _cigar.lengthSoftClippedLeft());
				_bases[pos].distFrom3 = TPseudoInt::fromLinear(f_m1_ps - pos);
			}
		}
	} else {
		// treat as single end
		if (_flags.isReverseStrand()) {
			// Single-end reverse
			//  3'  LL---x--RR  5' : readLength       = softclippedLeft + lengthSequenced + softClippedRight
			//  3'  ..---x--..  5' : lengthSequenced  = 6
			//  3'  LL...x....  5' : softClippledLeft = 2
			//  3'  .....x..RR  5' : softClippledRight= 2
			// pos: 0123456789
			// d5':   543210       : d5 = lengthSequenced - 1 - pos + softClippedLeft = 6 - 1 - 5 + 2 = 2
			// d3':   012345       : d3 = pos - softClippedLeft                       = 5 - 2         = 3
			for (size_t pos = _cigar.lengthSoftClippedLeft(); pos < length + _cigar.lengthSoftClippedLeft(); ++pos) {
				_bases[pos].distFrom5 = TPseudoInt::fromLinear(l_m1_ps - pos);
				_bases[pos].distFrom3 = TPseudoInt::fromLinear(pos - _cigar.lengthSoftClippedLeft());
			}
		} else {
			// Single-end forward
			//  5'  LL---x--RR  3' : readLength       = softclippedLeft + lengthSequenced + softClippedRight
			//  5'  ..---x--..  3' : lengthSequenced  = 6
			//  5'  LL...x....  3' : softClippledLeft = 2
			//  5'  .....x..RR  3' : softClippledRight= 2
			// pos: 0123456789
			// d5':   012345       : d5 = pos - softClippedLeft                       = 5 - 2         = 3
			// d3':   543210       : d3 = lengthSequenced - 1 - pos + softClippedLeft = 6 - 1 - 5 + 2 = 2
			for (size_t pos = _cigar.lengthSoftClippedLeft(); pos < length + _cigar.lengthSoftClippedLeft(); ++pos) {
				_bases[pos].distFrom5 = TPseudoInt::fromLinear(pos - _cigar.lengthSoftClippedLeft());
				_bases[pos].distFrom3 = TPseudoInt::fromLinear(l_m1_ps - pos);
			}
		}
	}
}

void TAlignment::_fillContext() {
	using genometools::Base;
	assert(!_bases.empty());
	if (_flags.isReverseStrand()) {
		// reverse
		for (size_t d = 0; d < _bases.size() - 1; ++d) {
			_bases[d].previousBase = _bases[d + 1].base;
		}
		_bases.back().previousBase = Base::N;
	} else {
		// forward
		_bases.front().previousBase = Base::N;
		for (size_t d = 1; d < _bases.size(); ++d)
			_bases[d].previousBase = _bases[d - 1].base;
	}
}

void TAlignment::parse() {
	_parseBasesQualities();
	_setQualitiesNoRecal();
}

void TAlignment::parse(const GenotypeLikelihoods::SequencingError::TModels &seqErrorModels) {
	_parseBasesQualities();

	// recalibrate
	seqErrorModels.recalibrate(*this);
	_sequenceAndQualitiesChanged = seqErrorModels.recalibrates();
}

void TAlignment::addReference(const genometools::TFastaReader &fasta) {
	const auto view = fasta.view(refID(), position(), _refSize);
	_referenceSequence.clear();
	std::copy(view.begin(), view.end(), std::back_inserter(_referenceSequence));
}

void TAlignment::setSequenceQualities(const TCigar &Cigar, const std::vector<genometools::Base> &Sequence,
									  const std::vector<PhredInt> &Qualities) {
	if (Cigar.lengthRead() != Sequence.size() || Cigar.lengthRead() != Qualities.size()) {
		DEVERROR("length of CIGAR, Sequences and Qualities do not match!");
	}
	_setCigar(Cigar);

	// parse bases and qualities
	_sequence.clear();
	_sequence.reserve(Sequence.size());
	_qualities.clear();
	_qualities.reserve(Sequence.size());
	for (size_t i = 0; i < Sequence.size(); ++i) {
		_sequence.push_back(genometools::base2char(Sequence[i]));
		_qualities.push_back(coretools::toChar(Qualities[i]));
	}
	_parseBasesQualities();
	_setQualitiesNoRecal();
	_sequenceAndQualitiesChanged = true; // will trigger that the strings are read form the bases
}

void TAlignment::setReadGroup(const uint16_t readGroupId) {
	_readGroupID = readGroupId;
}

//get overlap length from previous function and pass it on to cigar constructor
void TAlignment::merge(uint16_t overlapLength, size_t &mappedBasesClipped) {
	_setCigar(TCigar(cigar(), overlapLength, !isReverseStrand(), mappedBasesClipped));
}

//--------------------------------------
// getters
//--------------------------------------
bool TAlignment::isAlignedAtInternalPos(size_t internalPosition) const {
	assert(internalPosition < _alignedPosition.size());
	return _alignedPosition[internalPosition] >= 0;
}

uint32_t TAlignment::getLastInternalPos() const {
	return (_alignedPosition.size()-1);
}

genometools::Base TAlignment::referenceAtInternalPos(size_t internalPosition) const {
	assert(isAlignedAtInternalPos(internalPosition));
	assert((size_t)_alignedPosition[internalPosition] < _referenceSequence.size());
	return _referenceSequence[_alignedPosition[internalPosition]];
}

genometools::TGenomePosition TAlignment::positionInRef(size_t internalPosition) const {
	assert(isAlignedAtInternalPos(internalPosition));
	return *this + _alignedPosition[internalPosition];
}

void TAlignment::_updateSequenceAndQualities() const {
	if (_sequenceAndQualitiesChanged) {
		// update according to what is stored in bases
		_sequence.resize(_bases.size());
		_qualities.resize(_bases.size());

		for (size_t b = 0; b < _bases.size(); ++b) {
			_sequence[b]  = genometools::base2char(_bases[b].base);
			_qualities[b] = coretools::toChar(_bases[b].recalQuality);
		}

		_sequenceAndQualitiesChanged = false;
	}
}

std::string TAlignment::sequence() const {
	_updateSequenceAndQualities();
	return _sequence;
}

std::string TAlignment::qualities() const {
	_updateSequenceAndQualities();
	return _qualities;
}

//--------------------------------------------
// filters and other functions to modify data
//--------------------------------------------

void TAlignment::trimRead(uint64_t trimmingLength3Prime, uint64_t trimmingLength5Prime) {
	for (auto &b : _bases) {
		if (b.distFrom3.linear() < trimmingLength3Prime || b.distFrom5.linear() < trimmingLength5Prime) {
			b.base         = genometools::Base::N;
			b.recalQuality = PhredInt::highest();
		}
	}

	_sequenceAndQualitiesChanged = true;
}

void TAlignment::trimSoftClips() {
	// make sure read is parsed
	if (!_parsed) DEVERROR("Read was not parsed!");

	if (_cigar.lengthSoftClippedLeft()) {
		_bases.erase(_bases.begin(), _bases.begin() + _cigar.lengthSoftClippedLeft());
		_sequenceAndQualitiesChanged = true;
	}
	if (_cigar.lengthSoftClippedRight()) {
		_bases.resize(_bases.size() - _cigar.lengthSoftClippedRight());
		_sequenceAndQualitiesChanged = true;
	}
	_cigar.removeSoftClips();
}

void TAlignment::trimSoftClips(size_t maxNumberOfSoftClippedBases) {
	// make sure read is parsed
	if (!_parsed) DEVERROR("Read was not parsed!");

	// check if there is softclipping that exceeds the threshold on the left
	if (_cigar.lengthSoftClippedLeft() > maxNumberOfSoftClippedBases) {
		_bases.erase(_bases.begin(), _bases.begin() + (_cigar.lengthSoftClippedLeft() - maxNumberOfSoftClippedBases));
		_sequenceAndQualitiesChanged = true;
	}
	// then do the same on the right side
	if (_cigar.lengthSoftClippedRight() > maxNumberOfSoftClippedBases) {
		_bases.resize(_bases.size() - (_cigar.lengthSoftClippedRight() - maxNumberOfSoftClippedBases));
		_sequenceAndQualitiesChanged = true;
	}		
		// update cigar and length
	_cigar.removeSoftClips(maxNumberOfSoftClippedBases);
}


void TAlignment::binQualityScoresIllumina() {
	// make sure read is parsed
	if (!_parsed)
		DEVERROR("Read was not parsed!");

	// bin quality scores as done by Illumina
	for (auto &b : _bases) { b.recalQuality = makeIllumina(b.recalQuality); }

	_sequenceAndQualitiesChanged = true;
};

void TAlignment::recalibrateWithPMD(const GenotypeLikelihoods::TErrorModels &GLCalculator) {
	GLCalculator.recalibrateWithPMD(*this);
	_sequenceAndQualitiesChanged = true;
};

void TAlignment::setIsProperPair(bool ok) { _flags.setIsProperPair(ok); };

void TAlignment::downsampleAlignment(coretools::Probability fractionToKeep) {
	for (auto &b : _bases) {
		double r = coretools::instances::randomGenerator().getRand();
		if (r > fractionToKeep) {
			b.base                          = genometools::Base::N;
			b.recalQuality = PhredInt::highest();
		}
	}
	_sequenceAndQualitiesChanged = true;
};

}; // namespace BAM
