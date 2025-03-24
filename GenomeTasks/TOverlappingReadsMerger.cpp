#include "TOverlappingReadsMerger.h"
#include "TSequencedData.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include <memory>
#include <numeric>

namespace GenomeTasks {
using coretools::instances::logfile;

bool TOverlappingReadsMerger::_merge(BAM::TAlignment &Fwd, BAM::TAlignment &Rev) {
	if (Rev < Fwd) {
		//  FFFF
		// RRR
		++_cases[Cases::RStart_s_FStart];

		//  something strange, discard!
		return false;
	}

	const auto FStart = Fwd.position();
	const auto FLen   = Fwd.cigar().lengthMapped();
	const auto FEnd   = FStart + FLen;

	const auto RStart = Rev.position();
	const auto RLen   = Rev.cigar().lengthMapped();
	const auto REnd   = RStart + RLen;

	if (REnd < FEnd) {
		//  FFF
		// RRR
		++_cases[Cases::REnd_s_FEnd];

		//  something strange, discard!
		return false;
	}


	if (RStart >= FEnd) {
		// FFF    -> FFF
		//    RRR       RRR
		++_cases[Cases::NoOverlap];


		return true; // no overlap
	}

	// else:
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
	_merger->merge(Fwd, Rev, FEnd - RStart);
	return true;
}

void TOverlappingReadsMerger::_summary() {
	const double nCases = std::accumulate(_cases.begin(), _cases.end(), 0);
	logfile().startIndent("Merging summary:");
	logfile().list(_cases[Cases::NoOverlap], " pairs without overlap: ", _cases[Cases::NoOverlap]/nCases);
	logfile().list(_cases[Cases::Overlap], " pairs with overlap: ", _cases[Cases::Overlap]/nCases);
	logfile().list("Discarted ", _cases[Cases::BothFwd], " pairs with two Fwd strands: ", _cases[Cases::BothFwd]/nCases);
	logfile().list("Discarted ", _cases[Cases::BothRev], " pairs with two Rev strands: ", _cases[Cases::BothRev]/nCases);
	logfile().list("Discarted ", _cases[Cases::RStart_s_FStart], " pairs where Rev starts before Fwd: ", _cases[Cases::RStart_s_FStart]/nCases);
	logfile().list("Discarted ", _cases[Cases::REnd_s_FEnd], " pairs where Rev ends before Fwd: ", _cases[Cases::REnd_s_FEnd]/nCases);
	logfile().endIndent();
}

TOverlappingReadsMerger::TOverlappingReadsMerger() : TWaitingListBamTraverser("_merged.bam") {
	using coretools::instances::parameters;

	const auto method = parameters().get("mergingMethod", "middle");
	if (method == "random"){
		_merger = std::make_unique<TRandomMerger>();
		logfile().list("Merging method: ", method, ", will keep random read for all overlapping positions. (parameter 'mergingMethod')");
	}  else if(method == "keepFirst"){
		_merger = std::make_unique<TMateMerger>(BAM::Mate::first);
		logfile().list("Merging method: ", method, ", will keep read of first mate at overlapping positions. (parameter 'mergingMethod')");
	} else if(method == "keepFwd"){
		_merger = std::make_unique<TStrandMerger>(BAM::Strand::Fwd);
		logfile().list("Merging method: ", method, ", will keep read of forward strand at overlapping positions. (parameter 'mergingMethod')");
	} else if(method == "keepRev"){
		_merger = std::make_unique<TStrandMerger>(BAM::Strand::Rev);
		logfile().list("Merging method: ", method, ", will keep read of reversed strand at overlapping positions. (parameter 'mergingMethod')");
	} else if(method == "keepSecond"){
		_merger = std::make_unique<TMateMerger>(BAM::Mate::second);
		logfile().list("Merging method: ", method, ", will keep read of second mate at overlapping positions. (parameter 'mergingMethod')");
	} else if (method == "middle"){
		_merger = std::make_unique<TMiddleMerger>();
		logfile().list("Merging method: ", method, ", will keep half of the overlapping positions of each mate. (parameter 'mergingMethod')");
	} else {
		UERROR("Unknown merging method ", method, "! Use 'none', 'middle', 'keepFirst', 'keppSecond' or 'random'.");
	}

	if(!_genome.bamFile().filter(BAM::FilterType::ImproperPairs)){
		logfile().warning("Improper pairs are kept but will not be merged!");
	}
}

void TOverlappingReadsMerger::_handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) {
	// lhs <= rhs with respect to reference
	if (!lhs.alignment.isProperPair()) { // not a proper pair: mark mate as as improper too
		rhs.alignment.setIsProperPair(false);
		lhs.status = AlignmentStatus::filterOut;
		rhs.status = AlignmentStatus::filterOut;
	} else {
		// no parsing needed!
		auto merged = false;
		if (lhs.alignment.isReverseStrand() == rhs.alignment.isReverseStrand()) {
			if (lhs.alignment.isReverseStrand()) ++_cases[Cases::BothRev];
			else ++_cases[Cases::BothFwd];
			merged = false;
		} else {
			merged = rhs.alignment.isReverseStrand() ? _merge(lhs.alignment, rhs.alignment) : _merge(rhs.alignment, lhs.alignment);
		}

		if (merged) {
			lhs.status = AlignmentStatus::ready;
			rhs.status = AlignmentStatus::ready;
		} else {
			lhs.status = AlignmentStatus::filterOut;
			rhs.status = AlignmentStatus::filterOut;
		}

	}
}

void TOverlappingReadsMerger::run() {
	traverseBAM();
	_summary();
}

}
