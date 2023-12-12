#include "TBamFilter.h"
#include "TWaitingListBamTraverser.h"

namespace GenomeTasks {
using namespace coretools::str;

void TBamFilter::_handleMates(BAM::TAlignment &alignment, iterator mate) {
	if (!alignment.isProperPair()) {
		// not a proper pair: mark mate as as improper
		mate->alignment.setIsProperPair(false);
	}
	// mark both as ready for writing
	mate->status = AlignmentStatus::ready;
	_waitingList.emplace_back(alignment, AlignmentStatus::ready);
};

void TBamFilter::_handleSingle(BAM::TAlignment &alignment) {
	// read is single end: add for writing
	_waitingList.emplace_back(alignment, AlignmentStatus::ready);
};

bool TBamFilter::_alignmentCanBeWrittenUnchanged() {
	return !_recalibrate && !_genome.bamFile().curIsPaired() && _waitingList.empty() &&
		   (_removeSoftClippedBases
				? (_genome.bamFile().curCIGAR().lengthSoftClippedRight() < _maxNumberOfSoftClippedBases &&
				   _genome.bamFile().curCIGAR().lengthSoftClippedLeft() < _maxNumberOfSoftClippedBases)
				: true);
}

} // namespace GenomeTasks
