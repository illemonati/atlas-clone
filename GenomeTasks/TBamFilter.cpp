#include "TBamFilter.h"

namespace GenomeTasks {
using namespace coretools::str;

void TBamFilter::_handleMates(BAM::TAlignment &alignment, iterator mate) {
	if (!alignment.isProperPair()) {
		// not a proper pair: mark mate as as improper
		mate->setAsNonProperPair();
	}
	// mark both as ready for writing
	mate->makeReady();
	_alignmentStorage.emplace_back(&alignment, true);
};

void TBamFilter::_handleSingle(BAM::TAlignment &alignment) {
	// read is single end: add for writing
	_alignmentStorage.emplace_back(&alignment, true);
};

bool TBamFilter::_alignmentCanBeWrittenUnchanged() {
	return !_recalibrate && !_genome.bamFile().curIsPaired() && _alignmentStorage.empty() &&
		   (_removeSoftClippedBases
				? (_genome.bamFile().curCIGAR().lengthSoftClippedRight() < _maxNumberOfSoftClippedBases &&
				   _genome.bamFile().curCIGAR().lengthSoftClippedLeft() < _maxNumberOfSoftClippedBases)
				: true);
}

} // namespace GenomeTasks
