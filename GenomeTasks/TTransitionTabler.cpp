#include "TTransitionTabler.h"
#include "TSequencedData.h"
#include "coretools/Strings/toString.h"

namespace GenomeTasks {
using genometools::Base;
using coretools::str::toString;

TTransitionTabler::TTransitionTabler() {
	_parser.openReference(true);
	std::vector<std::string> header{"Strand", "Pos"};
	for (auto ref = Base::min; ref < Base::max; ++ref) {
		for (auto b = Base::min; b < Base::max; ++b) {
			header.push_back(coretools::str::toString(toString(ref), "->", toString(b)));
		}
	}
	_outTransitions[0].open(_genome.outputName() + "_transitionsSeq.txt.gz", header);
	_outTransitionsRel[0].open(_genome.outputName() + "_transitionsRelSeq.txt.gz", header);
	for (size_t m = 1; m < 3; ++m) {
		_outTransitions[m].open(toString(_genome.outputName(), "_transitions", m, ".txt.gz"), header);
		_outTransitionsRel[m].open(toString(_genome.outputName(), "_transitionsRel", m, ".txt.gz"), header);
	}
}

void TTransitionTabler::_handleAlignments(const BAM::TAlignment &first, const BAM::TAlignment &second) {
	if (first.isReverseStrand()) {
		_handleRev(_transitions[0][BAM::Strand::Rev], _transitions[1][BAM::Strand::Rev], first, 0);
	} else {
		_handleFwd(_transitions[0][BAM::Strand::Fwd], _transitions[1][BAM::Strand::Fwd], first, 0);
	}

	const auto len = first.cigar().lengthSequenced();
	if (second.isReverseStrand()) {
		_handleRev(_transitions[0][BAM::Strand::Rev], _transitions[2][BAM::Strand::Rev], second, len);
	} else {
		_handleFwd(_transitions[0][BAM::Strand::Fwd], _transitions[2][BAM::Strand::Fwd], second, len);
	}
}

void TTransitionTabler::_handleFwd(std::vector<Rho> &transAll, std::vector<Rho> &transMate, const BAM::TAlignment &aln,
								   size_t addL) {
	const auto cigar = aln.cigar();
	const auto len   = cigar.lengthAligned();
	const auto softL = cigar.lengthSoftClippedLeft();

	if (transAll.size() < len + addL) transAll.resize(len + addL);
	if (transMate.size() < len) transMate.resize(len);

	for (size_t i = 0; i < aln.size(); ++i) {
		if (aln.isAlignedAtInternalPos(i) && aln[i].base != Base::N && !aln[i].get<BAM::Flags::SoftClipped>()) {
			assert(i >= softL);
			assert(i - softL + addL < transAll.size());
			assert(i - softL < transMate.size());
			++transAll[i - softL + addL][aln.referenceAtInternalPos(i)][aln[i].base];
			++transMate[i - softL][aln.referenceAtInternalPos(i)][aln[i].base];
		}
	}
}

void TTransitionTabler::_handleRev(std::vector<Rho> &transAll, std::vector<Rho> &transMate, const BAM::TAlignment &aln,
                                   size_t addL) {
	const auto &cigar  = aln.cigar();
	const auto len     = cigar.lengthAligned();
	const auto softL   = cigar.lengthSoftClippedLeft();
	const auto l_m1_ps = len - 1 + softL;

	if (transAll.size() < len + addL) transAll.resize(len + addL);
	if (transMate.size() < len) transMate.resize(len);
	for (size_t i = 0; i < aln.size(); ++i) {
		if (aln.isAlignedAtInternalPos(i) && aln[i].base != Base::N && !aln[i].get<BAM::Flags::SoftClipped>()) {
			assert(l_m1_ps >= i);
			assert(l_m1_ps - i + addL < transAll.size());
			assert(l_m1_ps - i < transMate.size());
			++transAll[l_m1_ps - i + addL][aln.referenceAtInternalPos(i)][aln[i].base];
			++transMate[l_m1_ps - i][aln.referenceAtInternalPos(i)][aln[i].base];
		}
	}
}

void TTransitionTabler::_handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) {
	if (!lhs.alignment.isProperPair()) return;

	if (!lhs.alignment.isParsed()) { lhs.alignment.parse(); }
	if (!rhs.alignment.isParsed()) { rhs.alignment.parse(); }

	lhs.alignment.addReference(_parser.reference());
	rhs.alignment.addReference(_parser.reference());

	if (rhs.alignment.isSecondMate()) {
		_handleAlignments(lhs.alignment, rhs.alignment);
	} else {
		_handleAlignments(rhs.alignment, lhs.alignment);
	}
}
void TTransitionTabler::_handleSingle(TWaitingAlignment &) {
	// only interesting for paired-ended
}

void TTransitionTabler::run() {
	traverseBAM();
	_writeTransitions();
}

void TTransitionTabler::_writeTransitions() {
	using BAM::Strand;
	using genometools::Base;
	using genometools::base2char;
	for (auto strand = Strand::min; strand < Strand::max; ++strand) {
		coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base>, genometools::Base> rho{};
		for (size_t m = 0; m < 3; ++m) {
			auto &tr = _transitions[m][strand];
			for (size_t i = 0; i < tr.size(); ++i) {
				_outTransitions[m].write(strand, i);
				_outTransitionsRel[m].write(strand, i);
				coretools::TStrongArray<size_t, genometools::Base> tot{};
				for (auto ref = Base::min; ref < Base::max; ++ref) {
					// counts
					for (auto b = Base::min; b < Base::max; ++b) {
						_outTransitions[m].write(tr[i][ref][b]);
						tot[ref] += tr[i][ref][b];
						rho[ref][b] += tr[i][ref][b];
					}
					for (auto b = Base::min; b < Base::max; ++b) {
						// Rel
						const auto rel = tot[ref] ? double(tr[i][ref][b]) / tot[ref] : 0;
						_outTransitionsRel[m].write(fmt::format("{:.4f}", rel));
					}
				}
				_outTransitions[m].endln();
				_outTransitionsRel[m].endln();
			}
		}
	}
}
} // namespace GenomeTasks
