#ifndef TALIGNMENTTRAVERSER_H_
#define TALIGNMENTTRAVERSER_H_

#include "TAlignment.h"
#include "TParser.h"
#include "TReadTraverser.h"
#include "coretools/Main/TError.h"

namespace BAM {

class TAlignmentTraverser {
	TReadTraverser _readTraverser;
	TParser _parser;
	BAM::TAlignment _alignment;
	bool _filled = false;

	void _fill();

public:
	const BAM::TBamFile &bamFile() const noexcept { return _readTraverser.bamFile(); }
	BAM::TBamFile &bamFile() noexcept { return _readTraverser.bamFile(); }

	const BAM::RGInfo::TReadGroupInfo &rgInfo() const noexcept { return _readTraverser.rgInfo(); }
	BAM::RGInfo::TReadGroupInfo &rgInfo() noexcept { return _readTraverser.rgInfo(); }

	const GenotypeLikelihoods::TErrorModels &errorModels() const noexcept { return _readTraverser.errorModels(); };

	const std::string &outputName() const noexcept { return _readTraverser.outputName(); }

	const genometools::TChromosomes &chromosomes() const noexcept { return _readTraverser.chromosomes(); }
	const genometools::TChromosome &curChr() const noexcept { return _readTraverser.curChr(); }

	BAM::TAlignment &alignment();
	void jump(size_t RefID);

	void openReference(bool Required) { _parser.openReference(Required); }
	const genometools::TFastaReader& reference() const noexcept {return _parser.reference();};

	void nextAlignment();
	bool endOfAlignments() { return _readTraverser.endOfReads(); }
};
}

#endif
