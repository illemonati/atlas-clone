#include "TBamTraverser.h"

namespace GenomeTasks {
TFilteredBamTraverser::TFilteredBamTraverser() {
	const BAM::TBamFilters filters{true};
	_genome.bamFile().setFilters(filters);
}

void TFilteredBamTraverser::_traverseBAMPassedQC() {
	// parse through bam file
	_genome.bamFile().startProgressReporting();
	while (_genome.bamFile().readNextAlignmentThatPassesFilters()) {
		// handle alignment by derived classes
		_handleAlignment();

		// report
		_genome.bamFile().printProgress();
	}
	// report
	_genome.bamFile().printEndWithSummary(_genome.outputName());
}

TParsedBamTraverser::TParsedBamTraverser() {
	_genome.bamFile().setFilters(BAM::TBamFilters{true});
}


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
