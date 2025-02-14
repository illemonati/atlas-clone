/*
 * TAlignmentMerger.cpp
 *
 *  Created on: Oct 20, 2022
 *      Author: phaentu
 */

#include "TAlignmentMerger.h"

#include "TAlignment.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"

namespace GenomeTasks::AlignmentMerger {

using coretools::PhredInt;
using coretools::instances::randomGenerator;

namespace impl {
void callMergeFunction(BAM::TAlignment &alignment, BAM::TAlignment &mate, size_t overlapLength) {
	size_t lengthMappedBeforeMerge = alignment.cigar().lengthMapped();
	size_t mappedBasesClipped      = 0;
	// the reads are merged by softclipping bases, for which a new cigar string is constructed
	alignment.merge(overlapLength, mappedBasesClipped);
	if (alignment.isReverseStrand()) {
		// if the reverse strand gets merged, its position in the BAM-file as well as the position of the mate of the
		// first read need to be adjusted to account for the length of the added softclips on left side
		if (overlapLength < lengthMappedBeforeMerge) {
			alignment.moveOnRef(alignment.position() + mappedBasesClipped);
			mate.setMateGenomicPosition(alignment);
		} else {
			// if the overlap is larger than the read being merged, the position can only be moved up by the amount of
			// aligned bases in the read
			alignment.moveOnRef(alignment.position() + lengthMappedBeforeMerge);
			mate.setMateGenomicPosition(alignment);
		}
	}
}

std::pair<PhredInt, PhredInt> minQual(const BAM::TAlignment &firstRead, const BAM::TAlignment &secondRead) {
	// base iterator starts at first position of the reverse strand, then increments until it reaches either the last
	// aligned position of itself or the forward read
	std::vector<BAM::TSequencedData>::const_iterator dataIterator = secondRead.begin();
	size_t internalPos                                            = 0;
	// use the recalibrated quality
	PhredInt secondReadMinQual                                    = dataIterator->recalQuality;

	while (!secondRead.isAlignedAtInternalPos(internalPos) ||
	       (secondRead.positionInRef(internalPos).position() !=
	            firstRead.lastAlignedPositionWithRespectToRef().position() &&
	        secondRead.positionInRef(internalPos).position() !=
	            secondRead.lastAlignedPositionWithRespectToRef().position())) {
		// if the base at the current position of the iterator is aligned and its PhredIntProbability is higher
		// (therefore the error-probability is higher and the quality is lower) save this number as the new minimum
		// quality
		if (secondRead.isAlignedAtInternalPos(internalPos)) {
			if (dataIterator->recalQuality > secondReadMinQual) secondReadMinQual = dataIterator->recalQuality;
		}
		dataIterator++;
		internalPos++;
	}
	// base iterator starts at last position of the forward strand, then decrements until it reaches either the first
	// aligned position of itself or the forward read
	auto dataIteratorReverse  = firstRead.rbegin();
	internalPos               = firstRead.getLastInternalPos();
	PhredInt firstReadMinQual = dataIteratorReverse->recalQuality;
	while (!firstRead.isAlignedAtInternalPos(internalPos) ||
	       (firstRead.positionInRef(internalPos).position() != secondRead.position() &&
	        firstRead.positionInRef(internalPos).position() != firstRead.position())) {
		if (firstRead.isAlignedAtInternalPos(internalPos)) {
			if (dataIteratorReverse->recalQuality > firstReadMinQual)
				firstReadMinQual = dataIteratorReverse->recalQuality;
		}
		dataIteratorReverse++;
		internalPos--;
	}
	// returns a pair of the minimum qualities of both reads (in the overlap)
	return std::make_pair(firstReadMinQual, secondReadMinQual);
}

std::pair<PhredInt, PhredInt> getMinQuals(const BAM::TAlignment &alignment, const BAM::TAlignment &mate) {
	// this if-statement ensures that the order of the pair of minimum qualities that is returned matches the order in
	// which the two reads were passed to the function this is necessary because the minQual function needs the two
	// reads in the order of (forwardStrand, reverseStrand)
	if (!alignment.isReverseStrand()) {
		return minQual(alignment, mate);
	} else {
		std::pair<PhredInt, PhredInt> flippedResult = minQual(mate, alignment);
		return std::make_pair(flippedResult.second, flippedResult.first);
	}
}

PhredInt determineQualAtSingleBase(const BAM::TAlignment &alignment, size_t numBasesFromEnd, bool isFirst) {
	// for the first read, we start at the last base of the read and go to the base at the center of the overlap and
	// return its quality-value
	if (isFirst) {
		std::vector<BAM::TSequencedData>::const_reverse_iterator dataIterator = alignment.rbegin() + numBasesFromEnd;
		return dataIterator->recalQuality;
	} else {
		// for the second read, we start at the first base of the read
		std::vector<BAM::TSequencedData>::const_iterator dataIterator = alignment.begin() + numBasesFromEnd;
		return dataIterator->recalQuality;
	}
}

void compareQualities(const PhredInt &alignmentQual, const PhredInt &mateQual, size_t &alignmentOverlapLength,
                      size_t &mateOverlapLength) {
	if (alignmentQual < mateQual)
		mateOverlapLength++;
	else if (alignmentQual > mateQual)
		alignmentOverlapLength++;
	else {
		// if both qualities are the same, pick one randomly
		if (randomGenerator().pickOneOfTwo())
			mateOverlapLength++;
		else
			alignmentOverlapLength++;
	}
}

void checkIfReverseStrandFirst(const BAM::TAlignment &alignment, const BAM::TAlignment &mate,
                               size_t &alignmentOverlapLength, size_t &mateOverlapLength, bool alignmentIsFirst) {
	// if the reverse strand is the first read
	// we need to add the distance between the start positions to the overlap length of the reverse strand
	// and the distance between the end positions to the overlap length of the forward strand
	if (alignmentIsFirst && alignment.isReverseStrand()) {
		alignmentOverlapLength += mate.position() - alignment.position();
		mateOverlapLength +=
		    (mate.position() + mate.cigar().lengthMapped()) - (alignment.position() + alignment.cigar().lengthMapped());
	} else if (!alignmentIsFirst && mate.isReverseStrand()) {
		mateOverlapLength += alignment.position() - mate.position();
		alignmentOverlapLength +=
		    (alignment.position() + alignment.cigar().lengthMapped()) - (mate.position() + mate.cigar().lengthMapped());
	}
}

void mergeOddOverlap(BAM::TAlignment &firstRead, BAM::TAlignment &secondRead, size_t halfOverlap) {
	// calculate qualities for the position at the center of the overlap for both reads
	PhredInt firstReadQual  = determineQualAtSingleBase(firstRead, halfOverlap, true);
	PhredInt secondReadQual = determineQualAtSingleBase(secondRead, halfOverlap, false);

	size_t firstReadOverlapLength  = halfOverlap;
	size_t secondReadOverlapLength = halfOverlap;
	// compare qualities and add 1 to the overlap length of the read with the lower quality base at the center of the
	// overlap
	compareQualities(firstReadQual, secondReadQual, firstReadOverlapLength, secondReadOverlapLength);

	// check if the reverse strand comes first and adjust the overlap lengths accordingly if that is the case
	checkIfReverseStrandFirst(firstRead, secondRead, firstReadOverlapLength, secondReadOverlapLength, true);
	// merge
	callMergeFunction(firstRead, secondRead, firstReadOverlapLength);
	callMergeFunction(secondRead, firstRead, secondReadOverlapLength);
}

void mergeOverlapLargerThanRead(BAM::TAlignment &largerRead, BAM::TAlignment &smallerRead,
                                size_t largerReadOverlapLength, size_t smallerReadOverlapLength) {
	// if the overlap is divisible by two, you can just merge
	if (smallerRead.cigar().lengthMapped() % 2 == 0) {
		callMergeFunction(smallerRead, largerRead, smallerReadOverlapLength);
		callMergeFunction(largerRead, smallerRead, largerReadOverlapLength);
	} else {
		// if the overlap is not divisible by two, we first need to determine the quality at the position in the center
		// of the overlap for both reads
		PhredInt smallerReadQual =
		    determineQualAtSingleBase(smallerRead, smallerReadOverlapLength, !smallerRead.isReverseStrand());
		PhredInt largerReadQual =
		    determineQualAtSingleBase(largerRead, largerReadOverlapLength, !largerRead.isReverseStrand());

		compareQualities(smallerReadQual, largerReadQual, smallerReadOverlapLength, largerReadOverlapLength);
		callMergeFunction(smallerRead, largerRead, smallerReadOverlapLength);
		callMergeFunction(largerRead, smallerRead, largerReadOverlapLength);
	}
}

void handleOverlapLargerThanRead(BAM::TAlignment &largerRead, BAM::TAlignment &smallerRead) {
	// since the larger read fully eclipses the smaller read, the overlap is now as large as the number of aligned
	// positions in the smaller read
	size_t halfOverlap              = smallerRead.cigar().lengthMapped() / 2;
	size_t smallerReadOverlapLength = halfOverlap;
	size_t largerReadOverlapLength  = halfOverlap;
	// if the larger read is the forward strand, add the distance between the ends of both reads to its overlap
	// if it is the reverse read, add the distance between both starts
	if (!largerRead.isReverseStrand()) {
		largerReadOverlapLength += (largerRead.position() + largerRead.cigar().lengthMapped()) -
		                           (smallerRead.position() + smallerRead.cigar().lengthMapped());
	} else {
		largerReadOverlapLength += smallerRead.position() - largerRead.position();
	}

	// if both reads are forward/reverse, the smaller read needs to be set to reverse/forward so different sides of the
	// reads are soft-clipped
	if ((largerRead.isReverseStrand() && smallerRead.isReverseStrand())) {
		smallerRead.setIsReverseStrand(false);
		mergeOverlapLargerThanRead(largerRead, smallerRead, largerReadOverlapLength, smallerReadOverlapLength);
		smallerRead.setIsReverseStrand(true);
	} else if ((!largerRead.isReverseStrand() && !smallerRead.isReverseStrand())) {
		smallerRead.setIsReverseStrand(true);
		mergeOverlapLargerThanRead(largerRead, smallerRead, largerReadOverlapLength, smallerReadOverlapLength);
		smallerRead.setIsReverseStrand(false);
	} else {
		mergeOverlapLargerThanRead(largerRead, smallerRead, largerReadOverlapLength, smallerReadOverlapLength);
	}
}
} // namespace impl

//-----------------------------------------
// TAlignmentMergerType
//-----------------------------------------
size_t TBase::_merge(BAM::TAlignment &alignment, BAM::TAlignment &mate) {
	bool alignmentFirst  = mate > alignment;
	size_t overlapLength = 0;
	if ((alignment.isReverseStrand() && mate.isReverseStrand())) {
		alignmentFirst ? alignment.setIsReverseStrand(false) : mate.setIsReverseStrand(false);
		overlapLength = overlapLengthAndMerge(alignment, mate);
		!alignment.isReverseStrand() ? alignment.setIsReverseStrand(true) : mate.setIsReverseStrand(true);
	} else if ((!alignment.isReverseStrand() && !mate.isReverseStrand())) {
		alignmentFirst ? mate.setIsReverseStrand(true) : alignment.setIsReverseStrand(true);
		overlapLength = overlapLengthAndMerge(alignment, mate);
		alignment.isReverseStrand() ? alignment.setIsReverseStrand(false) : mate.setIsReverseStrand(false);
	} else {
		overlapLength = overlapLengthAndMerge(alignment, mate);
	}
	return overlapLength;
};

size_t TBase::determineOverlapLength(const BAM::TAlignment &alignment, const BAM::TAlignment &mate) {
	// No overlap if either one of the reads is already completely softclipped
	if (alignment.cigar().lengthAligned() == 0 || mate.cigar().lengthAligned() == 0) return 0;
	// if-statements check for overlap using the first aligned position in the BAM-file and the aligned length given in
	// the cigar string
	if (alignment.isReverseStrand()) {
		if (alignment.position() < mate.position() + mate.cigar().lengthMapped()) {
			size_t overlapLength = mate.position() + mate.cigar().lengthMapped() - alignment.position();
			return overlapLength;
		}
		return 0;
	} else {
		if (mate.position() < alignment.position() + alignment.cigar().lengthMapped()) {
			size_t overlapLength = alignment.position() + alignment.cigar().lengthMapped() - mate.position();
			return overlapLength;
		}
		return 0;
	}
}

size_t TBase::overlapLengthAndMerge(BAM::TAlignment &alignment, BAM::TAlignment &mate) {
	// check if reads overlap
	size_t overlapLength = determineOverlapLength(alignment, mate);
	// if they do -> merge
	if (overlapLength > 0) { impl::callMergeFunction(alignment, mate, overlapLength); }
	return overlapLength;
}

// TAlignmentMergerType_randomRead
//---------------------------------
TRandomRead::TRandomRead() : TBase() {};

size_t TRandomRead::merge(BAM::TAlignment &alignment, BAM::TAlignment &mate) {
	// randomly either merge mate or alignment
	if (randomGenerator().pickOneOfTwo()) return _merge(mate, alignment);
	// else
	return _merge(alignment, mate);
};

// TAlignmentMergerType_middle
//---------------------------------
TMiddle::TMiddle() : TBase() {};

size_t TMiddle::merge(BAM::TAlignment &alignment, BAM::TAlignment &mate) {
	const std::pair<size_t, bool> overlapLength = determineOverlapLength(alignment, mate);
	if (overlapLength.first > 0) {
		// in case one of the reads fully covers its mate
		if (overlapLength.first >= alignment.cigar().lengthMapped()) {
			impl::handleOverlapLargerThanRead(mate, alignment);
			return overlapLength.first;
		} else if (overlapLength.first >= mate.cigar().lengthMapped()) {
			impl::handleOverlapLargerThanRead(alignment, mate);
			return overlapLength.first;
		}
		// if both reads are forward/reverse strands, we need to change one of their directions before merging so we
		// don't clip the same side twice
		if ((alignment.isReverseStrand() && mate.isReverseStrand())) {
			!overlapLength.second ? mate.setIsReverseStrand(false) : alignment.setIsReverseStrand(false);
			sameDirectionMerge(alignment, mate, overlapLength);
			!alignment.isReverseStrand() ? alignment.setIsReverseStrand(true) : mate.setIsReverseStrand(true);
		} else if ((!alignment.isReverseStrand() && !mate.isReverseStrand())) {
			!overlapLength.second ? alignment.setIsReverseStrand(true) : mate.setIsReverseStrand(true);
			sameDirectionMerge(alignment, mate, overlapLength);
			alignment.isReverseStrand() ? alignment.setIsReverseStrand(false) : mate.setIsReverseStrand(false);
		} else {
			size_t halfOverlap = overlapLength.first / 2;
			// if the overlap length is divisible by 2, we only need to check if the reverse strand is first before
			// merging
			if (overlapLength.first % 2 == 0) {
				size_t alignmentOverlapLength = halfOverlap;
				size_t mateOverlapLength      = halfOverlap;
				impl::checkIfReverseStrandFirst(alignment, mate, alignmentOverlapLength, mateOverlapLength,
				                                overlapLength.second);

				impl::callMergeFunction(alignment, mate, alignmentOverlapLength);
				impl::callMergeFunction(mate, alignment, mateOverlapLength);
			} else {
				// this if-statement is necessary because the mergeOddOverlap-function needs the alignments in the order
				// of (firstRead, secondRead) to calculate the qualities at the center position
				if (overlapLength.second == true)
					impl::mergeOddOverlap(alignment, mate, halfOverlap);
				else
					impl::mergeOddOverlap(mate, alignment, halfOverlap);
			}
		}
	}
	return overlapLength.first;
};

std::pair<size_t, bool> TMiddle::determineOverlapLength(const BAM::TAlignment &alignment,
                                                                        const BAM::TAlignment &mate) {
	// this function works like the one it overrides, except that it checks which read comes first instead of which is
	// the reverse strand also it returns a boolean which states whether the alignment is the first read
	if (alignment.cigar().lengthAligned() == 0 || mate.cigar().lengthAligned() == 0) return std::make_pair(0, true);
	if (alignment > mate) {
		if (alignment.position() < mate.position() + mate.cigar().lengthMapped()) {
			size_t overlapLength = mate.position() + mate.cigar().lengthMapped() - alignment.position();
			bool alignmentFirst  = false;
			return std::make_pair(overlapLength, alignmentFirst);
		}
		// else
		return std::make_pair(0, false);
	} else if (mate > alignment) {
		if (mate.position() < alignment.position() + alignment.cigar().lengthMapped()) {
			size_t overlapLength = alignment.position() + alignment.cigar().lengthMapped() - mate.position();
			bool alignmentFirst  = true;
			return std::make_pair(overlapLength, alignmentFirst);
		}
		// else
		return std::make_pair(0, true);
	} else {
		// if both reads have the same start position, the overlap is the lower mapped length of the two
		size_t overlapLength;
		(alignment.cigar().lengthMapped() > mate.cigar().lengthMapped())
		    ? overlapLength = mate.cigar().lengthMapped()
		    : overlapLength = alignment.cigar().lengthMapped();
		bool alignmentFirst = randomGenerator().pickOneOfTwo();
		return std::make_pair(overlapLength, alignmentFirst);
	}
}

void TMiddle::sameDirectionMerge(BAM::TAlignment &alignment, BAM::TAlignment &mate,
                                                 std::pair<size_t, bool> overlapLength) {
	size_t halfOverlap = overlapLength.first / 2;
	if (overlapLength.first % 2 == 0) {
		size_t alignmentOverlapLength = halfOverlap;
		size_t mateOverlapLength      = halfOverlap;

		impl::callMergeFunction(alignment, mate, alignmentOverlapLength);
		impl::callMergeFunction(mate, alignment, mateOverlapLength);
	} else {
		// this if-statement is necessary because the mergeOddOverlap-function needs the alignments in the order of
		// (firstRead, secondRead) to calculate the qualities at the center position
		if (overlapLength.second == true)
			impl::mergeOddOverlap(alignment, mate, halfOverlap);
		else
			impl::mergeOddOverlap(mate, alignment, halfOverlap);
	}
}

// TAlignmentMergerType_firstMate
//---------------------------------
TFirstMate::TFirstMate() : TBase() {};

size_t TFirstMate::merge(BAM::TAlignment &alignment, BAM::TAlignment &mate) {
	// always merge first mate, keep second mate
	if (alignment.mate() == BAM::Mate::second)
		return _merge(mate, alignment);
	else
		return _merge(alignment, mate);
};

// TAlignmentMergerType_secondMate
//---------------------------------
TSecondMate::TSecondMate() : TBase() {};

size_t TSecondMate::merge(BAM::TAlignment &alignment, BAM::TAlignment &mate) {
	// always merge second mate, keep first mate
	if (alignment.mate() == BAM::Mate::second)
		return _merge(alignment, mate);
	else
		return _merge(mate, alignment);
};

// TAlignmentMergerType_highestQuality
//---------------------------------
THighestQuality::THighestQuality() : TBase() {};

size_t THighestQuality::merge(BAM::TAlignment &alignment, BAM::TAlignment &mate) {
	bool alignmentFirst  = mate > alignment;
	size_t overlapLength = 0;
	if ((alignment.isReverseStrand() && mate.isReverseStrand())) {
		alignmentFirst ? alignment.setIsReverseStrand(false) : mate.setIsReverseStrand(false);
		overlapLength = _overlapLengthAndMergeBetterQuality(alignment, mate);
		!alignment.isReverseStrand() ? alignment.setIsReverseStrand(true) : mate.setIsReverseStrand(true);
	} else if ((!alignment.isReverseStrand() && !mate.isReverseStrand())) {
		alignmentFirst ? mate.setIsReverseStrand(true) : alignment.setIsReverseStrand(true);
		overlapLength = _overlapLengthAndMergeBetterQuality(alignment, mate);
		alignment.isReverseStrand() ? alignment.setIsReverseStrand(false) : mate.setIsReverseStrand(false);
	} else {
		overlapLength = _overlapLengthAndMergeBetterQuality(alignment, mate);
	}
	return overlapLength;
}

size_t THighestQuality::_overlapLengthAndMergeBetterQuality(BAM::TAlignment &alignment, BAM::TAlignment &mate) {
	// check if reads overlap
	size_t overlapLength = TBase::determineOverlapLength(alignment, mate);
	// if they do -> calculate minimum quality of each read in the overlap
	if (overlapLength > 0) {
		std::pair<PhredInt, PhredInt> minQuals = impl::getMinQuals(alignment, mate);
		// merge the read with the lower minimum quality (the read with the higher PhredIntProbability has the higher
		// error-rate and the lower quality)
		if (minQuals.first < minQuals.second) {
			impl::callMergeFunction(mate, alignment, overlapLength);
		} else if (minQuals.second < minQuals.first) {
			impl::callMergeFunction(alignment, mate, overlapLength);
		} else {
			// if both reads have equal minimum qualities, randomly choose one
			if (randomGenerator().pickOneOfTwo()) {
				impl::callMergeFunction(mate, alignment, overlapLength);
			} else {
				impl::callMergeFunction(alignment, mate, overlapLength);
			}
		}
	}
	return overlapLength;
}

} // namespace GenomeTasks::AlignmentMerger
