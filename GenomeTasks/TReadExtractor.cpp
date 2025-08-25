#include "TReadExtractor.h"

#include "TBamFile.h"
#include "coretools/Main/TParameters.h"
#include "genometools/GenomePositions/TGenomeWindow.h"
#include "genometools/TAlleles.h"

namespace GenomeTasks {
using coretools::instances::parameters;
using coretools::instances::logfile;

void TReadExtractor::_traverseAlignments() {
	for (; !_alnTraverser.endOfAlignments(); _alnTraverser.nextAlignment()) {
		if (_alIt == _alleles.end()) { break; }

		auto &Alignment = _alnTraverser.alignment();

		// go forward
		while (Alignment.from() > _alIt->position) {
			++_alIt;
			if (_alIt == _alleles.end()) break;
		}
		DEBUG_ASSERT(_alIt->position.refID() >= _alnTraverser.alignment().refID());

		if (_alIt->position.refID() > _alnTraverser.alignment().refID()) _alnTraverser.jump(_alIt->position.refID());
		if (_alIt->position > Alignment.to()) continue;

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
		} else if (nRef && nAlt) {
			_nBoth++;
		} else if (nRef) {
			_nRef++;
			_outBAM.writeAlignment(Alignment);
		} else if (nAlt) {
			_nAlt++;
			Alignment.setReadGroup(Alignment.readGroupId() + _alnTraverser.bamFile().readGroups().size());
			_outBAM.writeAlignment(Alignment);
		}
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
	: _outBAM(_alnTraverser.outputName() + "_extr.bam", _alnTraverser.bamFile().samHeader(), _alnTraverser.bamFile().chromosomes(),
			  impl::mkRGs(_alnTraverser.bamFile().readGroups())) {}

void TReadExtractor::run() {
	_alleles.parse(parameters().get("alleles"), _alnTraverser.bamFile().chromosomes());
	_alIt = _alleles.begin();

	_traverseAlignments();
	logfile().list("N_Ref: ", _nRef, ", N_Alt: ", _nAlt, ", N_both: ", _nBoth, ", N_other: ", _nOther);
}
} // namespace GenomeTasks
