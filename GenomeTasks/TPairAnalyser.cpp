#include "TPairAnalyser.h"
#include "TSequencedData.h"

namespace GenomeTasks {

namespace impl {
using BAM::Mate;
using BAM::Strand;
constexpr Mate mate(bool IsSecondMate) {
	constexpr std::array<Mate, 2> mates = {Mate::first, Mate::second};
	return mates[IsSecondMate];
}

constexpr Strand strand(bool IsSecondMate) {
	constexpr std::array<Strand, 2> strands = {Strand::Fwd, Strand::Rev};
	return strands[IsSecondMate];
}
}

TPairAnalyser::TPairAnalyser() {
	_out.open(_genome.outputName() + "_pairs.txt.gz", {"Name", "Chr", "PosStart", "PosEnd", "Mate", "Strand", "cigar"});
}

void TPairAnalyser::_handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) {
	_out.writeNoDelim("\n"); // trick to prevent column-counter
	const auto& chrs = _genome.bamFile().chromosomes();
	for (const auto &aln : {lhs.alignment, rhs.alignment}) {
		_genome.bamFile().chromosomes()[lhs.alignment.refID()];
		_out.writeln(aln.name(), chrs[aln.refID()].name(), aln.position(), aln.position() + aln.cigar().lengthMapped(),
					 impl::mate(aln.isSecondMate()), impl::strand(aln.isReverseStrand()), aln.cigar().compileString());
	}
}
} // namespace GenomeTasks
