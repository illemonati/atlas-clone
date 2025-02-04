#include "TTransitionTabler.h"
#include "TSequencedData.h"

namespace GenomeTasks {
using genometools::Base;

TTransitionTabler::TTransitionTabler() {
	_parser.openReference(true);
	std::vector<std::string> header{"Strand", "Pos"};
	for (auto ref = Base::min; ref < Base::max; ++ref) {
		for (auto b = Base::min; b < Base::max; ++b) {
			header.push_back(coretools::str::toString(toString(ref), "->", toString(b)));
		}
	}
	_outTransitions.open(_genome.outputName() + "_transitions.txt.gz", header);
	_outTransitionsRel.open(_genome.outputName() + "_transitionsRel.txt.gz", header);
	_outTransitionsPsi.open(_genome.outputName() + "_transitionsPsi.txt.gz", header);
	_outTransitionsRho.open(_genome.outputName() + "_transitionsRho.txt.gz", header);
}
void TTransitionTabler::_handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) {
	if (lhs.alignment.isProperPair()) {
		if (!lhs.alignment.isParsed()) { lhs.alignment.parse(); }
		if (!rhs.alignment.isParsed()) { rhs.alignment.parse(); }
		lhs.alignment.addReference(_parser.reference());
		rhs.alignment.addReference(_parser.reference());

		const auto &first  = rhs.alignment.isSecondMate() ? lhs.alignment : rhs.alignment;
		const auto &second = rhs.alignment.isSecondMate() ? rhs.alignment : lhs.alignment;

		if (first.isReverseStrand()) {
			auto &trans1      = _transitions[BAM::Strand::Rev];
			const auto l_m1_ps = first.cigar().lengthSequenced() - 1 + first.cigar().lengthSoftClippedLeft();
			if (trans1.size() < first.parsedLength()) trans1.resize(first.parsedLength());
			for (size_t i = 0; i < first.parsedLength(); ++i) {
				if (first.isAlignedAtInternalPos(i) && first[i].base != Base::N &&
					!first[i].get<BAM::Flags::SoftClipped>()) {
					assert(l_m1_ps >= i);
					++trans1[l_m1_ps - i][first.referenceAtInternalPos(i)][first[i].base];
				}
			}
		} else {
			auto &trans1      = _transitions[BAM::Strand::Fwd];
			const auto softL1 = first.cigar().lengthSoftClippedLeft();
			if (trans1.size() < first.parsedLength()) trans1.resize(first.parsedLength());
			for (size_t i = 0; i < first.parsedLength(); ++i) {
				if (first.isAlignedAtInternalPos(i) && first[i].base != Base::N &&
					!first[i].get<BAM::Flags::SoftClipped>()) {
					assert(i >= softL1);
					++trans1[i - softL1][first.referenceAtInternalPos(i)][first[i].base];
				}
			}
		}
		const auto length1 = first.cigar().lengthSequenced();
		if (second.isReverseStrand()) {
			auto &trans2      = _transitions[BAM::Strand::Rev];
			const auto l_m1_ps = second.cigar().lengthSequenced() - 1 + second.cigar().lengthSoftClippedLeft();
			if (trans2.size() < second.parsedLength() + length1) trans2.resize(second.parsedLength() + length1);
			for (size_t i = 0; i < second.parsedLength(); ++i) {
				if (second.isAlignedAtInternalPos(i) && second[i].base != Base::N &&
					!second[i].get<BAM::Flags::SoftClipped>()) {
					assert(l_m1_ps >= i);
					++trans2[l_m1_ps - i + length1][second.referenceAtInternalPos(i)][second[i].base];
				}
			}
		} else {
			auto &trans2       = _transitions[BAM::Strand::Fwd];
			const auto softL2  = second.cigar().lengthSoftClippedLeft();
			if (trans2.size() < second.size() + length1) trans2.resize(first.size() + length1);
			for (size_t i = 0; i < second.parsedLength(); ++i) {
				if (second.isAlignedAtInternalPos(i) && second[i].base != Base::N &&
				    !second[i].get<BAM::Flags::SoftClipped>()) {
					assert(i >= softL2);
					++trans2[i - softL2 + length1][second.referenceAtInternalPos(i)][second[i].base];
				}
			}
		}
	}
}
void TTransitionTabler::_handleSingle(TWaitingAlignment &lhs) {
		if (!lhs.alignment.isParsed()) { lhs.alignment.parse(); }
		lhs.alignment.addReference(_parser.reference());

		const auto &single  = lhs.alignment;
		if (single.isReverseStrand()) {
			auto &trans      = _transitions[BAM::Strand::Rev];
			const auto l_m1_ps = single.cigar().lengthSequenced() - 1 + single.cigar().lengthSoftClippedLeft();
			if (trans.size() < single.parsedLength()) trans.resize(single.parsedLength());
			for (size_t i = 0; i < single.parsedLength(); ++i) {
				if (single.isAlignedAtInternalPos(i) && single[i].base != Base::N &&
					!single[i].get<BAM::Flags::SoftClipped>()) {
					assert(l_m1_ps >= i);
					++trans[l_m1_ps - i][single.referenceAtInternalPos(i)][single[i].base];
				}
			}
		} else {
			auto &trans      = _transitions[BAM::Strand::Fwd];
			const auto softL1 = single.cigar().lengthSoftClippedLeft();
			if (trans.size() < single.parsedLength()) trans.resize(single.parsedLength());
			for (size_t i = 0; i < single.parsedLength(); ++i) {
				if (single.isAlignedAtInternalPos(i) && single[i].base != Base::N &&
					!single[i].get<BAM::Flags::SoftClipped>()) {
					assert(i >= softL1);
					++trans[i - softL1][single.referenceAtInternalPos(i)][single[i].base];
				}
			}
		}
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
		auto &tr = _transitions[strand];
		for (size_t i = 0; i < tr.size(); ++i) {
			_outTransitions.write(strand, i);
			_outTransitionsRel.write(strand, i);
			_outTransitionsPsi.write(strand, i);
			_outTransitionsRho.write(strand, i);
			coretools::TStrongArray<size_t, genometools::Base> tot{};
			for (auto ref = Base::min; ref < Base::max; ++ref) {
				// counts
				for (auto b = Base::min; b < Base::max; ++b) {
					_outTransitions.write(tr[i][ref][b]);
					tot[ref] += tr[i][ref][b];
					rho[ref][b] += tr[i][ref][b];
				}
				const auto totRho = tot[ref] - tr[i][ref][ref];

				for (auto b = Base::min; b < Base::max; ++b) {
					// Rel
					const auto rel = tot[ref] ? double(tr[i][ref][b]) / tot[ref] : 0;
					_outTransitionsRel.write(fmt::format("{:.4f}", rel));

					// Rho
					if (ref == b) {
						_outTransitionsRho.write("  -  ");
					} else {
						const auto rho = totRho ? double(tr[i][ref][b]) / totRho : 0;
						_outTransitionsRho.write(fmt::format("{:.4f}", rho));
					}
				}
			}
			// Psi
			for (auto ref = Base::min; ref < Base::max; ++ref) {
				for (auto b = Base::min; b < Base::max; ++b) {
					if (ref == b) {
						_outTransitionsPsi.write("0.000");
					} else {
						const auto fromTo = tot[ref] ? double(tr[i][ref][b]) / tot[ref] : 0;
						const auto toFrom = tot[b] ? double(tr[i][b][ref]) / tot[b] : 0;
						const auto Psi    = (fromTo - toFrom) / (1.0 - toFrom);
						_outTransitionsPsi.write(fmt::format("{:.4f}", Psi));
					}
				}
			}
			_outTransitions.endln();
			_outTransitionsRel.endln();
			_outTransitionsRho.endln();
			_outTransitionsPsi.endln();
		}
	}

	// Flush once per chromosome
	_outTransitions.flush();
	_outTransitionsRel.flush();
	_outTransitionsRho.flush();
	_outTransitionsPsi.flush();
}
} // namespace GenomeTasks
