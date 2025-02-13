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
	_outTransitionsRho[0].open(_genome.outputName() + "_transitionsRhoSeq.txt.gz", header);
	for (size_t m = 1; m < 3; ++m) {
		_outTransitions[m].open(toString(_genome.outputName(), "_transitions", m, ".txt.gz"), header);
		_outTransitionsRel[m].open(toString(_genome.outputName(), "_transitionsRel", m, ".txt.gz"), header);
		_outTransitionsRho[m].open(toString(_genome.outputName() + "_transitionsRho", m, ".txt.gz"), header);
	}
}
void TTransitionTabler::_handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) {
	if (lhs.alignment.isProperPair()) {
		if (!lhs.alignment.isParsed()) { lhs.alignment.parse(); }
		if (!rhs.alignment.isParsed()) { rhs.alignment.parse(); }
		lhs.alignment.addReference(_parser.reference());
		rhs.alignment.addReference(_parser.reference());

		const auto &first  = rhs.alignment.isSecondMate() ? lhs.alignment : rhs.alignment;
		const auto &second = rhs.alignment.isSecondMate() ? rhs.alignment : lhs.alignment;
		const auto &cigar1 = first.cigar();
		const auto &cigar2 = second.cigar();
		const auto length1 = cigar1.lengthSequenced();
		const auto length2 = cigar2.lengthSequenced();

		if (first.isReverseStrand()) {
			auto &trans1  = _transitions[0][BAM::Strand::Rev];
			auto &transM1 = _transitions[1][BAM::Strand::Rev];

			const auto l_m1_ps = length1 - 1 + cigar1.lengthSoftClippedLeft();

			if (trans1.size() < length1) trans1.resize(length1);
			if (transM1.size() < length1) transM1.resize(length1);
			for (size_t i = 0; i < first.size(); ++i) {
				if (first.isAlignedAtInternalPos(i) && first[i].base != Base::N &&
					!first[i].get<BAM::Flags::SoftClipped>()) {
					assert(l_m1_ps >= i);
					++trans1[l_m1_ps - i][first.referenceAtInternalPos(i)][first[i].base];
					++transM1[l_m1_ps - i][first.referenceAtInternalPos(i)][first[i].base];
				}
			}
		} else {
			auto &trans1      = _transitions[0][BAM::Strand::Fwd];
			auto &transM1     = _transitions[1][BAM::Strand::Fwd];
			const auto softL1 = cigar1.lengthSoftClippedLeft();

			if (trans1.size() < length1) trans1.resize(length1);
			if (transM1.size() < length1) transM1.resize(length1);
			for (size_t i = 0; i < first.size(); ++i) {
				if (first.isAlignedAtInternalPos(i) && first[i].base != Base::N &&
					!first[i].get<BAM::Flags::SoftClipped>()) {
					assert(i >= softL1);
					++trans1[i - softL1][first.referenceAtInternalPos(i)][first[i].base];
					++transM1[i - softL1][first.referenceAtInternalPos(i)][first[i].base];
				}
			}
		}
		if (second.isReverseStrand()) {
			auto &trans2       = _transitions[0][BAM::Strand::Rev];
			auto &transM2      = _transitions[2][BAM::Strand::Rev];
			const auto l_m1_ps = length2 - 1 + cigar2.lengthSoftClippedLeft();

			if (trans2.size() < length2 + length1) trans2.resize(length2 + length1);
			if (transM2.size() < length2) transM2.resize(length2);
			for (size_t i = 0; i < second.size(); ++i) {
				if (second.isAlignedAtInternalPos(i) && second[i].base != Base::N &&
					!second[i].get<BAM::Flags::SoftClipped>()) {
					assert(l_m1_ps >= i);
					++trans2[l_m1_ps - i + length1][second.referenceAtInternalPos(i)][second[i].base];
					++transM2[l_m1_ps - i][second.referenceAtInternalPos(i)][second[i].base];
				}
			}
		} else {
			auto &trans2      = _transitions[0][BAM::Strand::Fwd];
			auto &transM2     = _transitions[2][BAM::Strand::Fwd];
			const auto softL2 = cigar2.lengthSoftClippedLeft();

			if (trans2.size() < length2 + length1) trans2.resize(length2 + length1);
			if (transM2.size() < length2 + length1) transM2.resize(length2);

			for (size_t i = 0; i < second.size(); ++i) {
				if (second.isAlignedAtInternalPos(i) && second[i].base != Base::N &&
				    !second[i].get<BAM::Flags::SoftClipped>()) {
					assert(i >= softL2);
					++trans2[i - softL2 + length1][second.referenceAtInternalPos(i)][second[i].base];
					++transM2[i - softL2][second.referenceAtInternalPos(i)][second[i].base];
				}
			}
		}
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
				_outTransitionsRho[m].write(strand, i);
				coretools::TStrongArray<size_t, genometools::Base> tot{};
				for (auto ref = Base::min; ref < Base::max; ++ref) {
					// counts
					for (auto b = Base::min; b < Base::max; ++b) {
						_outTransitions[m].write(tr[i][ref][b]);
						tot[ref] += tr[i][ref][b];
						rho[ref][b] += tr[i][ref][b];
					}
					const auto totRho = tot[ref] - tr[i][ref][ref];

					for (auto b = Base::min; b < Base::max; ++b) {
						// Rel
						const auto rel = tot[ref] ? double(tr[i][ref][b]) / tot[ref] : 0;
						_outTransitionsRel[m].write(fmt::format("{:.4f}", rel));

						// Rho
						if (ref == b) {
							_outTransitionsRho[m].write("  -  ");
						} else {
							const auto rho = totRho ? double(tr[i][ref][b]) / totRho : 0;
							_outTransitionsRho[m].write(fmt::format("{:.4f}", rho));
						}
					}
				}
				_outTransitions[m].endln();
				_outTransitionsRel[m].endln();
				_outTransitionsRho[m].endln();
			}
		}
	}
}
} // namespace GenomeTasks
