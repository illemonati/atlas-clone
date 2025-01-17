#include "TOverlappingReadsMerger.h"

#include "TAlignmentMerger.h"
#include "coretools/Main/TParameters.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;

TOverlappingReadsMerger::TOverlappingReadsMerger() : TWaitingListBamTraverser("_merged.bam") {
	if(!_genome.bamFile().filter(BAM::FilterType::ImproperPairs)){
		logfile().warning("Improper pairs are kept but will not be merged!");
	}

	//set merging method
	std::string method = parameters().get<std::string>("mergingMethod", "middle");
	if (method == "randomRead"){
		_merger = std::make_unique<AlignmentMerger::TRandomRead>();
		logfile().list("Merging method: will keep random read for all overlapping positions. (parameter 'mergingMethod')");
	} else if(method == "highestQuality"){
		_merger = std::make_unique<AlignmentMerger::THighestQuality>();
		logfile().list("Merging method: will keep read with highest minimum quality at overlapping positions. (parameter 'mergingMethod')");
	} else if(method == "firstMate"){
		_merger = std::make_unique<AlignmentMerger::TFirstMate>();
		logfile().list("Merging method: will merge first mate at overlapping positions. (parameter 'mergingMethod')");
	} else if(method == "secondMate"){
		_merger = std::make_unique<AlignmentMerger::TSecondMate>();
		logfile().list("Merging method: will merge second mate at overlapping positions. (parameter 'mergingMethod')");
	} else if (method == "middle"){
		_merger = std::make_unique<AlignmentMerger::TMiddle>();
		logfile().list("Merging method: will keep half of the overlapping positions of each mate. (parameter 'mergingMethod')");
	} else {
		UERROR("Unknown merging method ", method, "! Use 'none', 'middle', 'firstMate', 'secondMate', 'randomRead' or 'highestQuality'.");
	}
};

void TOverlappingReadsMerger::_handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) {
	if (!lhs.alignment.isProperPair()) { // not a proper pair: mark mate as as improper too
		rhs.alignment.setIsProperPair(false);
		lhs.status = AlignmentStatus::orphan;
		rhs.status = AlignmentStatus::orphan;
	} else {
		// attempt merging: make sure alignments are parsed
		// Note: if we recalibrate, they were already parsed
		if (!lhs.alignment.isParsed()) { lhs.alignment.parse(); }
		if (!rhs.alignment.isParsed()) { rhs.alignment.parse(); }

		_merger->merge(lhs.alignment, rhs.alignment);
		lhs.status = AlignmentStatus::ready;
		rhs.status = AlignmentStatus::ready;
	}
}

void TOverlappingReadsMerger::_handleSingle(TWaitingAlignment &lhs) {
	// add as ready for writing
	lhs.status = AlignmentStatus::ready;
};

}
