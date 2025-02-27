#include "TNewMerger.h"
#include "genometools/GenomePositions/TGenomePosition.h"

namespace GenomeTasks {

void TNewMerger::_mergeMiddle(BAM::TAlignment& left, BAM::TAlignment& right) {
	if (left.refID() != right.refID()) return;

	const auto lStart = left.position();
	const auto lLen   = left.cigar().lengthAligned();
	const auto lEnd   = lStart + lLen;

	const auto rStart = right.position();
	const auto rLen   = right.cigar().lengthAligned();
	const auto rEnd   = rStart + rLen;

	if (rStart >= lEnd) {
		// lll    -> lll
		//    rrr       rrr

		return; // no overlap
	}

	if (rStart > lStart && rEnd <= lEnd) {
		// llll -> llll
		//  rr      ss 

		right.cigar().setAllToSoftClipped();
		right.moveOnRef(-1);

	} else {
		//  lll  -> lls
		//   rrr     srr
		// or
		//  lll    -> lls
		//  rrrr      ssrrr

		const auto overlap  = lEnd - rStart;
		const auto lOverlap = overlap / 2;
		const auto rOverlap = overlap - lOverlap; // this takes care of odd numbers

		left.cigar().addSoftClipsRight(lOverlap);
		right.cigar().addSoftClipsLeft(rOverlap);

		right += rOverlap;
		left.setMateGenomicPosition(right);
	}
}

void TNewMerger::_merge(BAM::TAlignment &lhs, BAM::TAlignment &rhs) {
	// middleMerge
	if (lhs < rhs) {
		_mergeMiddle(lhs, rhs);
	} else if (rhs < lhs) {
		_mergeMiddle(rhs, lhs);
	} else {
		// same starting position -> shorter first
		if (lhs.cigar().lengthAligned() < rhs.cigar().lengthAligned()) {
			_mergeMiddle(lhs, rhs);
		} else {
			_mergeMiddle(rhs, lhs);
		}
	}
}

void TNewMerger::_handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) {
	if (!lhs.alignment.isProperPair()) { // not a proper pair: mark mate as as improper too
		rhs.alignment.setIsProperPair(false);
		lhs.status = AlignmentStatus::orphan;
		rhs.status = AlignmentStatus::orphan;
	} else {
		// no parsing needed!
		_merge(lhs.alignment, rhs.alignment);


		lhs.status = AlignmentStatus::ready;
		rhs.status = AlignmentStatus::ready;
	}
}

}
