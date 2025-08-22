#include "TBamTraverser.h"

namespace GenomeTasks {

void TParsedBamTraverser::_traverseBAMPassedQC() {
	// parse through bam file
	_genome.bamFile().startProgressReporting();
	BAM::TAlignment alignment;
	while (_genome.bamFile().readNextAlignmentThatPassesFilters()) {
		_parser.fill(_genome, alignment);
		_handleAlignment(alignment);

		// report
		_genome.bamFile().printProgress();
	}

	// report
	_genome.bamFile().printEndWithSummary(_genome.outputName());
}
}
