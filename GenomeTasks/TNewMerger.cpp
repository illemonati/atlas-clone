#include "TNewMerger.h"
#include "coretools/Main/TLog.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include <memory>
#include <numeric>

namespace GenomeTasks {

bool TMiddleMerger::_mergeMiddle(BAM::TAlignment &Fwd, BAM::TAlignment &Rev) {
	const auto FStart = Fwd.position();
	const auto FLen   = Fwd.cigar().lengthMapped();
	const auto FEnd   = FStart + FLen;

	const auto RStart = Rev.position();
	const auto RLen   = Rev.cigar().lengthMapped();
	const auto REnd   = RStart + RLen;

	if (RStart >= FEnd) {
		// FFF    -> FFF
		//    RRR       RRR
		++_cases[Cases::NoOverlap];

		return true; // no overlap
	}

	if (RStart < FStart) {
		//  FFFF
		// RRR
		++_cases[Cases::RStart_s_FStart];

		//  something strange, discard!
		return false;
	}

	if (REnd < FEnd) {
		//  FFF
		// RRR
		++_cases[Cases::REnd_s_FEnd];

		//  something strange, discard!
		return false;
	}

	//  FFF  -> FFs
	//   RRR     sRR
	// or
	//  FFFF -> FFFs
	//   RRR     ssR
	// or
	//  FFF    -> FFs
	//  RRRR      ssRR
	// or
	//  FFFF   -> FFss
	//  RRRR      ssRR
	++_cases[Cases::Overlap];

	const auto overlap  = FEnd - RStart;
	const auto FOverlap = overlap / 2;
	const auto ROverlap = overlap - FOverlap; // this takes care of odd numbers

	Fwd.cigar().addSoftClipsRight(FOverlap);
	Rev.cigar().addSoftClipsLeft(ROverlap);

	Rev += ROverlap;
	Fwd.setMateGenomicPosition(Rev);

	return true;
}

bool TMiddleMerger::merge(BAM::TAlignment &lhs, BAM::TAlignment &rhs) {
	if (lhs.isReverseStrand() == rhs.isReverseStrand()) {
		if (lhs.isReverseStrand()) ++_cases[Cases::BothRev];
		else ++_cases[Cases::BothFwd];

		return false;
	}

	// middleMerge
	if (rhs.isSecondMate()) {
		return _mergeMiddle(lhs, rhs);
	} else {
		return _mergeMiddle(rhs, lhs);
	}
}

void TMiddleMerger::summary() {
	using coretools::instances::logfile;
	const double nCases = std::accumulate(_cases.begin(), _cases.end(), 0);
	logfile().startIndent("Merging summary:");
	logfile().list(_cases[Cases::NoOverlap], " pairs without overlap: ", _cases[Cases::NoOverlap]/nCases);
	logfile().list(_cases[Cases::Overlap], " pairs with overlap: ", _cases[Cases::Overlap]/nCases);
	logfile().list(_cases[Cases::BothFwd], " pairs with two Fwd strands: ", _cases[Cases::BothFwd]/nCases);
	logfile().list(_cases[Cases::BothRev], " pairs with two Rev strands: ", _cases[Cases::BothRev]/nCases);
	logfile().list(_cases[Cases::RStart_s_FStart], " pairs where Rev starts before Fwd: ", _cases[Cases::RStart_s_FStart]/nCases);
	logfile().list(_cases[Cases::REnd_s_FEnd], " pairs where Rev ends before Fwd: ", _cases[Cases::REnd_s_FEnd]/nCases);
	logfile().endIndent();
}

TNewMerger::TNewMerger() : TWaitingListBamTraverser("_merged.bam") {
	_merger = std::make_unique<TMiddleMerger>();
}

void TNewMerger::_handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) {
	// lhs <= rhs with respect to reference
	if (!lhs.alignment.isProperPair()) { // not a proper pair: mark mate as as improper too
		rhs.alignment.setIsProperPair(false);
		lhs.status = AlignmentStatus::orphan;
		rhs.status = AlignmentStatus::orphan;
	} else {
		// no parsing needed!
		bool merged = _merger->merge(lhs.alignment, rhs.alignment);

		if (merged) {
			lhs.status = AlignmentStatus::ready;
			rhs.status = AlignmentStatus::ready;
		} else {
			lhs.status = AlignmentStatus::filterOut;
			rhs.status = AlignmentStatus::filterOut;
		}

	}
}

void TNewMerger::run() {
	traverseBAM();
	_merger->summary();
}

}
