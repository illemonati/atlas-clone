#ifndef TALIGNMENTTRAVERSER_H_
#define TALIGNMENTTRAVERSER_H_

#include "TAlignment.h"
#include "TBamFile.h"
#include "TBamFilters.h"
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
	TAlignmentTraverser() = default;
	TAlignmentTraverser(std::string_view Name, size_t I) : _readTraverser(Name, true, I) {}
	TAlignmentTraverser(std::string_view Name, size_t I, const TBamFilters &Filters) : _readTraverser(Name, false, I) {
		bamFile().setFilters(Filters);
	}

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

	void setSilent() {_readTraverser.setSilent();}
};
}

#endif
