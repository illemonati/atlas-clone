#include "TBamFilter.h"
#include "TWaitingListBamTraverser.h"

namespace GenomeTasks {
using namespace coretools::str;

void TBamFilter::_handleMates(TWaitingAlignment &Mate) {
	if (!_waitingList.back().alignment.isProperPair()) {
		// not a proper pair: mark mate as as improper
		Mate.alignment.setIsProperPair(false);
	}
	// mark both as ready for writing
	Mate.status                = AlignmentStatus::ready;
	_waitingList.back().status = AlignmentStatus::ready;
};

void TBamFilter::_handleSingle() {
	// read is single end: add for writing
	_waitingList.back().status = AlignmentStatus::ready;
};

bool TBamFilter::_alignmentCanBeWrittenUnchanged() {
	return !_recalibrate && !_genome.bamFile().curIsPaired() && _waitingList.empty() &&
		   (_removeSoftClippedBases
				? (_genome.bamFile().curCIGAR().lengthSoftClippedRight() < _maxNumberOfSoftClippedBases &&
				   _genome.bamFile().curCIGAR().lengthSoftClippedLeft() < _maxNumberOfSoftClippedBases)
				: true);
}

} // namespace GenomeTasks
