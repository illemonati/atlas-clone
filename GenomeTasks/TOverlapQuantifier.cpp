#include "TOverlapQuantifier.h"

#include "TAlignmentMerger.h"

namespace GenomeTasks {

using coretools::instances::logfile;

TOverlapQuantifier::TOverlapQuantifier() {
	// filter bamFile
	_genome.bamFile().setFilters(BAM::TBamFilters{true});
}

void TOverlapQuantifier::run(){
	//prepare counter
	coretools::TCountDistributionVector overlapDist;
	//parse BAM file
	_genome.bamFile().startProgressReporting();
	while(_genome.bamFile().readNextAlignment()){
		//if on new chromosome, empty storage
		if(_genome.bamFile().chrChanged()){
			//clear storage
			_alignmentStorage.clear();
		}

		//check if first alignment in storage is too far away from current alignment
		//if yes, first alignment in storage is considered an orphan
		auto it = _alignmentStorage.begin();
		while(it != _alignmentStorage.end() && (_genome.bamFile().curPosition() - it->alignment) > _genome.bamFile().maxReadLength()){
			it = _alignmentStorage.erase(it);
		}

		//check if read passed filters and is proper pair
		if(_genome.bamFile().curPassedQC() && _genome.bamFile().curIsProperPair()){
			//parse alignment
			BAM::TAlignment alignment;
			_genome.bamFile().fill(alignment);

			//check if mate is in storage.
			auto mate = std::find_if(_alignmentStorage.begin(), _alignmentStorage.end(),
									 [alignment](const auto &wa) { return wa.alignment.name() == alignment.name(); });
			if (mate == _alignmentStorage.end()) {
				    // add alignment to storage and wait for mate
				_alignmentStorage.emplace_back(alignment, AlignmentStatus::waiting);
			} else {
				    // mate found
				    if (alignment.readGroupId() != mate->alignment.readGroupId()) {
					UERROR("Mates '", alignment.name(), "' are in different read groups!");
				}

				//calculate overlap and fragment length and add to storage
					size_t overlap = AlignmentMerger::TBase::determineOverlapLength(alignment, mate->alignment);
				size_t fragmentLength = alignment.fragmentLength();

				overlapDist.add(fragmentLength, overlap);
			}
		}

		//report
		_genome.bamFile().printProgress();
	}

	//done parsing bam file: report
	_genome.bamFile().printSummary(_genome.outputName());

	//write distribution
	std::string filename = _genome.outputName() + "_overlapStats.txt";
	logfile().listFlush("Writing distribution of fragment length and overlap to file '" + filename + "' ...");
	const std::vector<std::string> header = {"fragmentLength", "overlap", "count"};
	coretools::TOutputFile out(filename, header);
	overlapDist.write(out);
	out.close();
	logfile().done();
};

}
