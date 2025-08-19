#include "TReadExtractor.h"

#include "TBamFile.h"
#include "coretools/Main/TParameters.h"
#include "genometools/GenomePositions/TGenomeWindow.h"
#include "genometools/TAlleles.h"

namespace GenomeTasks {
using coretools::instances::parameters;
using coretools::instances::logfile;

void TReadExtractor::_handleAlignment(BAM::TAlignment& Alignment) {
	if (_alIt == _alleles.end()) {
		_genome.bamFile().jumpToEnd();
		return;
	}

	// go forward
	while (Alignment.from() > _alIt->position) {
		++_alIt;
		if (_alIt == _alleles.end()) return;
	}
	auto it = _alIt;

	size_t nRef   = 0;
	size_t nAlt   = 0;
	size_t nOther = 0;

	for (size_t i = 0; i < Alignment.size(); ++i) {
		if (!Alignment.isAlignedAtInternalPos(i)) continue;
		if (Alignment.positionInRef(i) == it->position) {
			const auto d = Alignment[i].base;
			if (d == it->ref) {
				++nRef;
			} else if (d == it->alt) {
				++nAlt;
			} else {
				++nOther;
			}
			++it;
			if (it == _alleles.end() || it->position > Alignment.to()) break;
		}
	}

	if (nOther) {
		_nOther++;
		return;
	}
	if (nRef && nAlt) {
		_nBoth ++;
		return;
	}

	if (nRef) {
		_nRef++;
		_outBAM.writeAlignment(Alignment);
	} else if (nAlt) {
		_nAlt++;
		Alignment.setReadGroup(Alignment.readGroupId() + _genome.bamFile().readGroups().size());
		_outBAM.writeAlignment(Alignment);
	}
}

namespace impl {
auto mkRGs(BAM::TReadGroups RGs) {
	const auto size = RGs.size();
	for (size_t i = 0; i < size; ++i) { RGs.addAlternativeRG(RGs[i].name_ID + "_alt", RGs[i].name_ID); }
	return RGs;
}
} // namespace impl

TReadExtractor::TReadExtractor()
	: _outBAM(_genome.outputName() + "_extr.bam", _genome.bamFile().samHeader(), _genome.bamFile().chromosomes(),
			  impl::mkRGs(_genome.bamFile().readGroups())) {}

void TReadExtractor::run() {
	_alleles.parse(parameters().get("alleles"), _genome.bamFile().chromosomes());
	_alIt = _alleles.begin();
	_genome.bamFile().jump(_alIt->position);

	_traverseBAMPassedQC();
	logfile().list("N_Ref: ", _nRef, ", N_Alt: ", _nAlt, ", N_both: ", _nBoth, ", N_other: ", _nOther);
}
} // namespace GenomeTasks
