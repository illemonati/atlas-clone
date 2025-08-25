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

	void _fill() {
		if (!_filled) {
			_parser.fill(_readTraverser, _alignment);
			_filled = true;
		}
	}

public:
	const BAM::TBamFile &bamFile() const noexcept { return _readTraverser.bamFile(); }
	BAM::TBamFile &bamFile() noexcept { return _readTraverser.bamFile(); }

	const BAM::RGInfo::TReadGroupInfo &rgInfo() const noexcept { return _readTraverser.rgInfo(); }
	BAM::RGInfo::TReadGroupInfo &rgInfo() noexcept { return _readTraverser.rgInfo(); }

	const GenotypeLikelihoods::TErrorModels &errorModels() const noexcept { return _readTraverser.errorModels(); };

	const std::string &outputName() const noexcept { return _readTraverser.outputName(); }

	const genometools::TChromosomes &chromosomes() const noexcept { return _readTraverser.chromosomes(); }
	const genometools::TChromosome &curChr() const noexcept { return _readTraverser.curChr(); }

	BAM::TAlignment &alignment() noexcept {
		_fill(); // only fill, which is time-consuming, if alignment is needed
		return _alignment;
	}

	void jump(size_t RefID) {
		bamFile().jump(RefID);
		nextAlignment();
	}

	void openReference(bool Required) {
		_parser.openReference(Required);
	}

	void nextAlignment() {
		DEBUG_ASSERT(!endOfAlignments());
		_readTraverser.nextRead();
		_filled = false;
	}

	bool endOfAlignments() {
		return _readTraverser.endOfReads();
	}
};	
}

#endif
