#include "TBamFilter.h"
#include "TWaitingListBamTraverser.h"

namespace GenomeTasks {
using namespace coretools::str;

void TBamFilter::_handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) {
	if (!lhs.alignment.isProperPair()) {
		// not a proper pair: mark mate as as improper
		rhs.alignment.setIsProperPair(false);
	}
	// mark both as ready for writing
	lhs.status = AlignmentStatus::ready;
	rhs.status = AlignmentStatus::ready;
}

void TBamFilter::_handleSingle(TWaitingAlignment &lhs) {
	// read is single end: add for writing
	lhs.status = AlignmentStatus::ready;
}

bool TBamFilter::_alignmentCanBeWrittenUnchanged() {
	return !_recalibrate && !_genome.bamFile().curIsPaired() && _waitingList.empty() &&
		   (_removeSoftClippedBases
				? (_genome.bamFile().curCIGAR().lengthSoftClippedRight() < _maxNumberOfSoftClippedBases &&
				   _genome.bamFile().curCIGAR().lengthSoftClippedLeft() < _maxNumberOfSoftClippedBases)
				: true);
}

} // namespace GenomeTasks
