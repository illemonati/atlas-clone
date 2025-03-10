#include "TMerger.h"
#include "coretools/Main/TRandomGenerator.h"
#include <type_traits>

namespace GenomeTasks {

namespace impl {
template<typename SomeInt>
constexpr bool isOdd(SomeInt N) {
	static_assert(std::is_integral_v<SomeInt>);
	return N & 1;
}
}

void TMiddleMerger::merge(BAM::TAlignment &Fwd, BAM::TAlignment &Rev, size_t Overlap) {
	// all edge cases are already handeld!

	// odd number -> one more for Mate1
	const auto FOverlap = impl::isOdd(Overlap) ? Overlap / 2 + !Fwd.isSecondMate() : Overlap / 2;
	const auto ROverlap = Overlap - FOverlap; // this takes care of odd numbers

	Fwd.cigar().addSoftClipsRight(FOverlap);
	Rev.cigar().addSoftClipsLeft(ROverlap);

	Rev.advanceOnRef(ROverlap);
	Fwd.setMateGenomicPosition(Rev.from());
}

void TRandomMerger::merge(BAM::TAlignment &Fwd, BAM::TAlignment &Rev, size_t Overlap) {
	using coretools::instances::randomGenerator;
	// all edge cases are already handeld!

	if (randomGenerator().getRand() < 0.5) {
		Fwd.cigar().addSoftClipsRight(Overlap);
	} else {
		Rev.cigar().addSoftClipsLeft(Overlap);

		Rev.advanceOnRef(Overlap);
		Fwd.setMateGenomicPosition(Rev.from());
	}
}

void TMateMerger::merge(BAM::TAlignment &Fwd, BAM::TAlignment &Rev, size_t Overlap) {
	// all edge cases are already handeld!

	if (Fwd.mate() != _keep) {
		Fwd.cigar().addSoftClipsRight(Overlap);
	} else {
		Rev.cigar().addSoftClipsLeft(Overlap);

		Rev.advanceOnRef(Overlap);
		Fwd.setMateGenomicPosition(Rev.from());
	}
}

void TStrandMerger::merge(BAM::TAlignment &Fwd, BAM::TAlignment &Rev, size_t Overlap) {
	// all edge cases are already handeld!

	if (_keep != BAM::Strand::Fwd) {
		Fwd.cigar().addSoftClipsRight(Overlap);
	} else {
		Rev.cigar().addSoftClipsLeft(Overlap);

		Rev.advanceOnRef(Overlap);
		Fwd.setMateGenomicPosition(Rev.from());
	}
}

}

